#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ==========================================
 * 物理 NAND 与系统级参数定义
 * ========================================== */
#define USER_DATA_BYTES 4096
#define USER_CRC_BYTES 4
#define META_DATA_BYTES 16
#define META_CRC_BYTES 4

#define LDPC_PAYLOAD_BYTES (USER_DATA_BYTES + USER_CRC_BYTES + META_DATA_BYTES + META_CRC_BYTES)
#define LDPC_PAYLOAD_BITS (LDPC_PAYLOAD_BYTES * 8) // 32960 bits

#define PARITY_BYTES 512
#define PARITY_BITS (PARITY_BYTES * 8) // 4096 bits
#define CODEWORD_BITS (LDPC_PAYLOAD_BITS + PARITY_BITS)
#define CODEWORD_BYTES (CODEWORD_BITS / 8) // 4632 Bytes

#define LIFTING_FACTOR 64
#define NUM_DATA_BLOCKS (LDPC_PAYLOAD_BITS / LIFTING_FACTOR)
#define NUM_CHK_BLOCKS (PARITY_BITS / LIFTING_FACTOR)
#define COLUMN_WEIGHT 4
#define MAX_ITERATIONS 20
#define MAX_HARD_ITERS 1000

#define SCRAMBLER_31_BIT 1

#define PAGES_PER_BLOCK 64
#define NUM_PHYSICAL_BLOCKS 32
#define NUM_PHYSICAL_PAGES (NUM_PHYSICAL_BLOCKS * PAGES_PER_BLOCK)
#define NUM_LOGICAL_PAGES 1024
#define INVALID_PPA 0xFFFFFFFF

typedef unsigned int LBA;
typedef unsigned int PPA;

/* ==========================================
 * 全局数据结构与物理模型
 * ========================================== */
struct ldpc_rom_entry
{
    unsigned short row;
    unsigned short col;
    unsigned short shift;
};
struct check_node
{
    int *connected_vns;
    int degree;
};
struct var_node
{
    int *connected_cns;
    int degree;
};

static struct ldpc_rom_entry *g_rom_base_matrix = NULL;
static int g_rom_matrix_entries = 0;
static struct check_node *g_check_nodes = NULL;
static struct var_node *g_var_nodes = NULL;

typedef struct
{
    float v_th;
} NandCell;

typedef struct
{
    int is_valid;
    int read_count;
    NandCell cells[CODEWORD_BITS];
} PhysicalPage;

typedef struct
{
    char inject_meta_mismatch;
    char inject_erase_error;
    char inject_unc_error;
    char inject_crc_error;
    char inject_all_zero_error;
} Inject_err_t;

static PhysicalPage g_nand_flash[NUM_PHYSICAL_PAGES];
static PPA g_l2p_table[NUM_LOGICAL_PAGES];
static int g_active_block = 0;
static int g_page_offset = 0;

/* ==========================================
 * 工具函数：Bit 与 Byte 转换
 * ========================================== */
void bytes_to_bits(const unsigned char *bytes, unsigned char *bits, int num_bytes)
{
    for (int i = 0; i < num_bytes; i++)
        for (int b = 0; b < 8; b++)
            bits[i * 8 + b] = (bytes[i] >> b) & 1;
}

void bits_to_bytes(const unsigned char *bits, unsigned char *bytes, int num_bytes)
{
    for (int i = 0; i < num_bytes; i++)
    {
        unsigned char val = 0;
        for (int b = 0; b < 8; b++)
            val |= (bits[i * 8 + b] << b);
        bytes[i] = val;
    }
}

unsigned int calculate_crc32(const unsigned char *data, int length)
{
    unsigned int crc = 0xFFFFFFFF;
    for (int i = 0; i < length; i++)
    {
        crc ^= data[i];
        for (int j = 0; j < 8; j++)
            crc = (crc >> 1) ^ (0xEDB88320 & (-(crc & 1)));
    }
    return ~crc;
}

void pack_nand_page(const unsigned char *user_data, const unsigned char *meta_data, unsigned char *ldpc_payload)
{
    int offset = 0;
    memcpy(ldpc_payload + offset, user_data, USER_DATA_BYTES);
    offset += USER_DATA_BYTES;
    unsigned int user_crc = calculate_crc32(user_data, USER_DATA_BYTES);
    memcpy(ldpc_payload + offset, &user_crc, USER_CRC_BYTES);
    offset += USER_CRC_BYTES;

    memcpy(ldpc_payload + offset, meta_data, META_DATA_BYTES);
    offset += META_DATA_BYTES;
    unsigned int meta_crc = calculate_crc32(meta_data, META_DATA_BYTES);
    memcpy(ldpc_payload + offset, &meta_crc, META_CRC_BYTES);
}

int verify_nand_page(const unsigned char *decoded_payload)
{
    int offset = 0;
    unsigned int calc_user_crc = calculate_crc32(decoded_payload + offset, USER_DATA_BYTES);
    unsigned int read_user_crc;
    memcpy(&read_user_crc, decoded_payload + offset + USER_DATA_BYTES, USER_CRC_BYTES);
    if (calc_user_crc != read_user_crc)
        return -1;
    offset += USER_DATA_BYTES + USER_CRC_BYTES;

    unsigned int calc_meta_crc = calculate_crc32(decoded_payload + offset, META_DATA_BYTES);
    unsigned int read_meta_crc;
    memcpy(&read_meta_crc, decoded_payload + offset + META_DATA_BYTES, META_CRC_BYTES);
    if (calc_meta_crc != read_meta_crc)
        return -2;

    return 0;
}

#if SCRAMBLER_31_BIT
void nand_data_randomizer(unsigned char *data, int length, unsigned int seed)
{
    unsigned int lfsr = (seed ^ 0x5AA55AA5) & 0x7FFFFFFF;
    if (lfsr == 0)
        lfsr = 1;
    for (int i = 0; i < length; i++)
    {
        unsigned char random_byte = 0;
        for (int b = 0; b < 8; b++)
        {
            unsigned int bit = ((lfsr >> 30) ^ (lfsr >> 27)) & 1;
            lfsr = ((lfsr << 1) | bit) & 0x7FFFFFFF;
            random_byte = (random_byte << 1) | bit;
        }
        data[i] ^= random_byte;
    }
}
#endif

/* ==========================================
 * 模块 1.8: 真实物理层读写引擎与可视化
 * ========================================== */
float generate_gaussian_noise(void)
{
    float u1 = ((float)rand() / RAND_MAX) + 1e-9f;
    float u2 = ((float)rand() / RAND_MAX) + 1e-9f;
    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
}

void physical_program_page(PPA ppa, const unsigned char *codeword_bits)
{
    float sigma = 0.4f;
    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        if (codeword_bits[i] == 0)
        {
            g_nand_flash[ppa].cells[i].v_th = 3.0f + generate_gaussian_noise() * sigma;
        }
        else
        {
            g_nand_flash[ppa].cells[i].v_th = -3.0f + generate_gaussian_noise() * sigma;
        }
    }
    g_nand_flash[ppa].read_count = 0;
    g_nand_flash[ppa].is_valid = 1;
}

// 可视化：打印物理 Vth 分布的 ASCII 柱状图
// 可视化：打印物理 Vth 分布的 ASCII 柱状图 (横坐标电压，纵坐标数量)
void print_vth_histogram(PPA ppa, const char *title)
{
    printf("\n=== Vth 分布柱状图 [%s] (PPA %04d) ===\n", title, ppa);

    const int NUM_BINS = 61; // 奇数，保证有一个绝对中心的 0V 刻度
    int bins[61] = {0};
    float min_v = -6.0f, max_v = 6.0f; // 恢复到标准观察区间
    float step = (max_v - min_v) / NUM_BINS;

    int count_0 = 0; // 记录物理电压 > 0V 的 Cell 数量 (Program态)
    int count_1 = 0; // 记录物理电压 < 0V 的 Cell 数量 (Erase态)

    // 1. 统计电压落入各个 Bin 的晶胞数量，并统计 0/1 比例
    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        float v = g_nand_flash[ppa].cells[i].v_th;

        // 以 0.0V 为理想的硬判决读取阈值 (V_read)
        if (v > 0.0f)
        {
            count_0++; // 导电沟道关闭，读出 0
        }
        else
        {
            count_1++; // 导电沟道开启，读出 1
        }

        int idx = (int)((v - min_v) / step);
        if (idx < 0)
            idx = 0;
        if (idx >= NUM_BINS)
            idx = NUM_BINS - 1;
        bins[idx]++;
    }

    // 2. 寻找最大峰值，用于动态拉伸 Y 轴高度
    int max_cnt = 0;
    for (int i = 0; i < NUM_BINS; i++)
    {
        if (bins[i] > max_cnt)
            max_cnt = bins[i];
    }

    // 3. 从上到下逐行打印图表
    const int CHART_HEIGHT = 15; // 控制图表的高度(行数)
    for (int h = CHART_HEIGHT; h > 0; h--)
    {
        int current_y = max_cnt * h / CHART_HEIGHT;
        printf("%4d | ", current_y); // 打印 Y 轴刻度

        for (int i = 0; i < NUM_BINS; i++)
        {
            int bar_h = (int)((float)bins[i] / max_cnt * CHART_HEIGHT);
            if (bar_h >= h)
            {
                float v_center = min_v + i * step;
                // 用 # 标出 0.0V 附近的极易混淆危险区
                if (v_center > -0.2f && v_center < 0.2f)
                {
                    printf("#");
                }
                else
                {
                    printf("*");
                }
            }
            else
            {
                printf(" ");
            }
        }
        printf("\n");
    }

    // 4. 打印底部的 X 轴横线
    printf("     +-");
    for (int i = 0; i < NUM_BINS; i++)
        printf("-");
    printf("\n");

    // 5. 打印 X 轴的电压刻度标签 (精确对齐)
    printf("Vth:  ");
    for (int i = 0; i < NUM_BINS; i++)
    {
        if (i == 0)
        {
            printf("-6V");
            i += 2;
        }
        else if (i == 15)
        {
            printf("-3V");
            i += 2;
        }
        else if (i == 30)
        {
            printf("0V");
            i += 1;
        }
        else if (i == 45)
        {
            printf("+3V");
            i += 2;
        }
        else if (i == 60)
        {
            printf("+6V");
            i += 2;
        }
        else
            printf(" ");
    }
    printf("\n");

    // 6. 【新增】：打印物理硬判决 0/1 比例统计
    float ratio_0 = (float)count_0 / CODEWORD_BITS * 100.0f;
    float ratio_1 = (float)count_1 / CODEWORD_BITS * 100.0f;
    printf("-----------------------------------------------------------------------\n");
    printf(" [物理硬判决统计 (V_read = 0.0V)]\n");
    printf(" -> 读出 '0' (Program 态, Vth > 0V): %5d bits (%5.2f%%)\n", count_0, ratio_0);
    printf(" -> 读出 '1' (Erase 态,   Vth < 0V): %5d bits (%5.2f%%)\n", count_1, ratio_1);
    printf("=======================================================================\n\n");
}

// 物理老化模拟：强制施加读干扰和热噪声
void force_age_page(PPA ppa, int read_cycles)
{
    g_nand_flash[ppa].read_count += read_cycles;
    float disturb_shift = g_nand_flash[ppa].read_count * 0.006f; // Read Disturb 导致擦除态往右漂移
    float noise_sigma = 0.6f;                                    // 老化后分布变宽

    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        float vth = g_nand_flash[ppa].cells[i].v_th;
        if (vth < 0)
            vth += disturb_shift;
        vth += generate_gaussian_noise() * noise_sigma;
        g_nand_flash[ppa].cells[i].v_th = vth;
    }
}

// 物理动作 4: 真实的双重编程 (Double Program) - 基于硅片实测数据修正
void force_double_program_page(PPA ppa)
{
    float sigma = 0.4f;

    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        float current_vth = g_nand_flash[ppa].cells[i].v_th;

        // 模拟向这个 Page 写入完全随机的新数据
        int new_bit = rand() % 2;

        if (current_vth < 0)
        {
            // 【情况 A】：原本是 Erase 态 (-3V)
            if (new_bit == 0)
            {
                // 一半被推到 Program 态 (+3V)
                g_nand_flash[ppa].cells[i].v_th = 3.0f + generate_gaussian_noise() * sigma;
            }
            else
            {
                // 另一半保持 Erase 态，叠加少量 Program Disturb 干扰
                g_nand_flash[ppa].cells[i].v_th += 0.3f + generate_gaussian_noise() * 0.1f;
            }
        }
        else
        {
            // 【情况 B】：原本已经是 Program 态 (+3V)
            // 根据实测：由于 ISPP Verify 保护和 FN 隧穿自限性，阈值电压基本不变
            // 无论是写 0 还是写 1，只会增加一点点随机的电荷扰动使其分布变宽
            g_nand_flash[ppa].cells[i].v_th += 0.1f + generate_gaussian_noise() * 0.15f;
        }
    }
}

const double g_llr_lut[8] = {15.0, 8.0, 0.0, 3.0, -15.0, -8.0, 0.0, -3.0};

void physical_soft_read_page(PPA ppa, double *llr_out, Inject_err_t err)
{
    g_nand_flash[ppa].read_count++;
    float read_noise_sigma = 0.1f;

    if (err.inject_unc_error)
    {
        read_noise_sigma = 2.0f;
    }

    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        float vth = g_nand_flash[ppa].cells[i].v_th;
        if (err.inject_erase_error)
            vth = -3.0f;
        if (err.inject_all_zero_error)
            vth = 3.0f;

        vth += generate_gaussian_noise() * read_noise_sigma;

        double rx_symbol = vth / 3.0;
        unsigned char hb = (rx_symbol < 0.0) ? 1 : 0;
        unsigned char sb1 = (rx_symbol > -0.2 && rx_symbol <= 0.2) ? 1 : 0;
        unsigned char sb2 = (rx_symbol > -0.6 && rx_symbol <= 0.6) ? 1 : 0;
        llr_out[i] = g_llr_lut[(hb << 2) | (sb1 << 1) | sb2];
    }
}

/* ==========================================
 * 模块 2: LDPC 引擎核心
 * ========================================== */
void build_mock_rom_table(void)
{
    int j, w, row_idx, shift_val, c_prime, r_prime, diff;
    int temp_matrix[NUM_CHK_BLOCKS][NUM_DATA_BLOCKS];
    int current_entry = 0, is_valid, attempts, restart_count = 0;

restart_matrix:
    if (restart_count > 5000)
        exit(-1);
    memset(temp_matrix, -1, sizeof(temp_matrix));
    current_entry = 0;

    for (j = 0; j < NUM_DATA_BLOCKS; j++)
    {
        for (w = 0; w < COLUMN_WEIGHT; w++)
        {
            int row_attempts = 0;
            do
            {
                row_idx = rand() % NUM_CHK_BLOCKS;
                if (++row_attempts > 100)
                    goto restart_matrix;
            } while (temp_matrix[row_idx][j] != -1);

            attempts = 0;
            do
            {
                shift_val = rand() % LIFTING_FACTOR;
                is_valid = 1;
                for (c_prime = 0; c_prime < j; c_prime++)
                {
                    if (temp_matrix[row_idx][c_prime] != -1)
                    {
                        for (r_prime = 0; r_prime < NUM_CHK_BLOCKS; r_prime++)
                        {
                            if (r_prime != row_idx && temp_matrix[r_prime][j] != -1 && temp_matrix[r_prime][c_prime] != -1)
                            {
                                diff = (temp_matrix[r_prime][c_prime] + shift_val) - (temp_matrix[r_prime][j] + temp_matrix[row_idx][c_prime]);
                                if ((diff % LIFTING_FACTOR + LIFTING_FACTOR) % LIFTING_FACTOR == 0)
                                {
                                    is_valid = 0;
                                    break;
                                }
                            }
                        }
                    }
                    if (!is_valid)
                        break;
                }
                if (++attempts > LIFTING_FACTOR * 2)
                {
                    restart_count++;
                    goto restart_matrix;
                }
            } while (!is_valid);
            temp_matrix[row_idx][j] = shift_val;
        }
    }

    g_rom_base_matrix = (struct ldpc_rom_entry *)malloc(NUM_DATA_BLOCKS * COLUMN_WEIGHT * sizeof(struct ldpc_rom_entry));
    for (j = 0; j < NUM_DATA_BLOCKS; j++)
    {
        for (int i = 0; i < NUM_CHK_BLOCKS; i++)
        {
            if (temp_matrix[i][j] != -1)
            {
                g_rom_base_matrix[current_entry].row = i;
                g_rom_base_matrix[current_entry].col = j;
                g_rom_base_matrix[current_entry].shift = temp_matrix[i][j];
                current_entry++;
            }
        }
    }
    g_rom_matrix_entries = current_entry;
}

int ldpc_engine_init(void)
{
    int i, entry_idx, z_offset, cn_idx, vn_idx;
    g_check_nodes = (struct check_node *)calloc(PARITY_BITS, sizeof(struct check_node));
    g_var_nodes = (struct var_node *)calloc(CODEWORD_BITS, sizeof(struct var_node));

    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            cn_idx = g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR + z_offset;
            vn_idx = g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR + (z_offset + g_rom_base_matrix[entry_idx].shift) % LIFTING_FACTOR;
            g_check_nodes[cn_idx].degree++;
            g_var_nodes[vn_idx].degree++;
        }
    }
    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].degree++;
        g_var_nodes[LDPC_PAYLOAD_BITS + i].degree++;
    }

    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].connected_vns = (int *)malloc(g_check_nodes[i].degree * sizeof(int));
        g_check_nodes[i].degree = 0;
    }
    for (i = 0; i < CODEWORD_BITS; i++)
    {
        g_var_nodes[i].connected_cns = (int *)malloc(g_var_nodes[i].degree * sizeof(int));
        g_var_nodes[i].degree = 0;
    }

    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            cn_idx = g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR + z_offset;
            vn_idx = g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR + (z_offset + g_rom_base_matrix[entry_idx].shift) % LIFTING_FACTOR;
            g_check_nodes[cn_idx].connected_vns[g_check_nodes[cn_idx].degree++] = vn_idx;
            g_var_nodes[vn_idx].connected_cns[g_var_nodes[vn_idx].degree++] = cn_idx;
        }
    }
    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].connected_vns[g_check_nodes[i].degree++] = LDPC_PAYLOAD_BITS + i;
        g_var_nodes[LDPC_PAYLOAD_BITS + i].connected_cns[g_var_nodes[LDPC_PAYLOAD_BITS + i].degree++] = i;
    }
    return 0;
}

void ldpc_encode(const unsigned char *data_in, unsigned char *codeword_out)
{
    int entry_idx, z_offset, data_bit_idx, parity_bit_idx;
    memcpy(codeword_out, data_in, LDPC_PAYLOAD_BITS);
    memset(codeword_out + LDPC_PAYLOAD_BITS, 0, PARITY_BITS);

    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            data_bit_idx = g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR + z_offset;
            parity_bit_idx = g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR + (z_offset + LIFTING_FACTOR - g_rom_base_matrix[entry_idx].shift) % LIFTING_FACTOR;
            codeword_out[LDPC_PAYLOAD_BITS + parity_bit_idx] ^= data_in[data_bit_idx];
        }
    }
}

static double extract_sign(double value) { return (value >= 0) ? 1.0 : -1.0; }

int ldpc_decode_soft(double *channel_llr, unsigned char *out_codeword)
{
    double **msg_chk_to_var = (double **)malloc(PARITY_BITS * sizeof(double *));
    int max_node_degree = 0, iter, i, j, edge_idx, search_idx, is_all_zero;
    for (i = 0; i < PARITY_BITS; i++)
    {
        msg_chk_to_var[i] = (double *)calloc(g_check_nodes[i].degree, sizeof(double));
        if (g_check_nodes[i].degree > max_node_degree)
            max_node_degree = g_check_nodes[i].degree;
    }
    double *posteriori_llr = (double *)malloc(CODEWORD_BITS * sizeof(double));
    double *msg_var_to_chk = (double *)malloc(max_node_degree * sizeof(double));

    for (iter = 0; iter < MAX_ITERATIONS; iter++)
    {
        for (i = 0; i < PARITY_BITS; i++)
        {
            double min_mag1 = 1e9, min_mag2 = 1e9, total_sign_prod = 1.0;
            int min_mag1_idx = -1;
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                int target_vn = g_check_nodes[i].connected_vns[edge_idx];
                msg_var_to_chk[edge_idx] = (iter == 0) ? channel_llr[target_vn] : (posteriori_llr[target_vn] - msg_chk_to_var[i][edge_idx]);
                double current_mag = fabs(msg_var_to_chk[edge_idx]);
                total_sign_prod *= extract_sign(msg_var_to_chk[edge_idx]);
                if (current_mag < min_mag1)
                {
                    min_mag2 = min_mag1;
                    min_mag1 = current_mag;
                    min_mag1_idx = edge_idx;
                }
                else if (current_mag < min_mag2)
                {
                    min_mag2 = current_mag;
                }
            }
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                double excl_sign = total_sign_prod * extract_sign(msg_var_to_chk[edge_idx]);
                double excl_mag = (edge_idx == min_mag1_idx) ? min_mag2 : min_mag1;
                msg_chk_to_var[i][edge_idx] = 0.75 * excl_sign * excl_mag;
            }
        }

        is_all_zero = 1;
        for (j = 0; j < CODEWORD_BITS; j++)
        {
            posteriori_llr[j] = channel_llr[j];
            for (edge_idx = 0; edge_idx < g_var_nodes[j].degree; edge_idx++)
            {
                int target_cn = g_var_nodes[j].connected_cns[edge_idx];
                int back_edge_idx = -1;
                for (search_idx = 0; search_idx < g_check_nodes[target_cn].degree; search_idx++)
                {
                    if (g_check_nodes[target_cn].connected_vns[search_idx] == j)
                    {
                        back_edge_idx = search_idx;
                        break;
                    }
                }
                posteriori_llr[j] += msg_chk_to_var[target_cn][back_edge_idx];
            }
            out_codeword[j] = (posteriori_llr[j] < 0) ? 1 : 0;
        }

        for (i = 0; i < PARITY_BITS; i++)
        {
            unsigned char current_syndrome = 0;
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
                current_syndrome ^= out_codeword[g_check_nodes[i].connected_vns[edge_idx]];
            if (current_syndrome != 0)
            {
                is_all_zero = 0;
                break;
            }
        }

        if (is_all_zero)
        {
            for (i = 0; i < PARITY_BITS; i++)
                free(msg_chk_to_var[i]);
            free(msg_chk_to_var);
            free(posteriori_llr);
            free(msg_var_to_chk);
            return iter + 1;
        }
    }
    for (i = 0; i < PARITY_BITS; i++)
        free(msg_chk_to_var[i]);
    free(msg_chk_to_var);
    free(posteriori_llr);
    free(msg_var_to_chk);
    return -1;
}

/* ==========================================
 * 模块 3: FTL 与全数据通路交互
 * ========================================== */
void ftl_init(void)
{
    for (int i = 0; i < NUM_LOGICAL_PAGES; i++)
        g_l2p_table[i] = INVALID_PPA;
    for (int i = 0; i < NUM_PHYSICAL_PAGES; i++)
    {
        g_nand_flash[i].is_valid = 0;
        g_nand_flash[i].read_count = 0;
        for (int j = 0; j < CODEWORD_BITS; j++)
        {
            g_nand_flash[i].cells[j].v_th = -3.0f; // 初始化为完美的物理擦除态
        }
    }
}

PPA allocate_physical_page(void)
{
    if (g_page_offset >= PAGES_PER_BLOCK)
    {
        g_active_block++;
        g_page_offset = 0;
    }
    if (g_active_block >= NUM_PHYSICAL_BLOCKS)
    {
        printf("[Fatal] 闪存已满！\n");
        exit(-1);
    }
    return g_active_block * PAGES_PER_BLOCK + g_page_offset++;
}

void ftl_write_full_path(LBA lba, const unsigned char *host_user_data)
{
    PPA new_ppa = allocate_physical_page();
    if (g_l2p_table[lba] != INVALID_PPA)
        g_nand_flash[g_l2p_table[lba]].is_valid = 0;

    unsigned char ftl_meta_data[META_DATA_BYTES] = {0};
    memcpy(ftl_meta_data, &lba, sizeof(LBA));

    unsigned char payload_bytes[LDPC_PAYLOAD_BYTES];
    pack_nand_page(host_user_data, ftl_meta_data, payload_bytes);

    nand_data_randomizer(payload_bytes, LDPC_PAYLOAD_BYTES, new_ppa);

    unsigned char payload_bits[LDPC_PAYLOAD_BITS];
    unsigned char codeword_bits[CODEWORD_BITS];
    bytes_to_bits(payload_bytes, payload_bits, LDPC_PAYLOAD_BYTES);

    ldpc_encode(payload_bits, codeword_bits);

    physical_program_page(new_ppa, codeword_bits);

    g_l2p_table[lba] = new_ppa;
    printf("[FW Write] LBA %04d 映射至 PPA %04d | 物理电荷注入完毕。\n", lba, new_ppa);
}

void ftl_read_full_path(LBA lba, Inject_err_t err)
{
    PPA target_ppa = g_l2p_table[lba];
    if (target_ppa == INVALID_PPA)
    {
        printf("[FW Read]  LBA %04d => [空数据]\n", lba);
        return;
    }

    if (err.inject_meta_mismatch)
        target_ppa = (target_ppa > 0) ? (target_ppa - 1) : (1);

    printf("[FW Read]  请求读取 LBA %04d (物理寻址: PPA %04d)...\n", lba, target_ppa);

    double rx_llr[CODEWORD_BITS];
    unsigned char rx_soft_out_bits[CODEWORD_BITS];

    physical_soft_read_page(target_ppa, rx_llr, err);

    int iter = ldpc_decode_soft(rx_llr, rx_soft_out_bits);
    if (iter < 0)
    {
        printf("  -> [Error] LDPC UE (不可纠错错误)，触发 RAIN 重构！\n");
        return;
    }

    unsigned char decoded_payload[LDPC_PAYLOAD_BYTES];
    bits_to_bytes(rx_soft_out_bits, decoded_payload, LDPC_PAYLOAD_BYTES);

    unsigned int descramble_seed = target_ppa;
    if (err.inject_crc_error)
        descramble_seed = target_ppa + 1;
    nand_data_randomizer(decoded_payload, LDPC_PAYLOAD_BYTES, descramble_seed);

    int crc_res = verify_nand_page(decoded_payload);
    if (crc_res == -1)
    {
        printf("  -> [Error] User CRC 校验失败 (静默数据损坏 SDC)！\n");
        return;
    }
    if (crc_res == -2)
    {
        printf("  -> [Error] Meta CRC 校验失败！\n");
        return;
    }

    LBA read_lba;
    memcpy(&read_lba, decoded_payload + USER_DATA_BYTES + USER_CRC_BYTES, sizeof(LBA));
    if (read_lba != lba)
    {
        printf("  -> [Error] 地址错乱！读出的 OOB 指向 LBA %d！\n", read_lba);
        return;
    }

    unsigned char final_host_data[16];
    memcpy(final_host_data, decoded_payload, 15);
    final_host_data[15] = '\0';

    printf("  -> [Success] 纠错耗时 %d 轮。提取数据: \"%s\"\n", iter, final_host_data);
}

void ftl_rw_test()
{
    printf("\n--- 第一阶段：主机初始化写入 ---\n");
    unsigned char host_buffer[USER_DATA_BYTES] = {0};

    strcpy((char *)host_buffer, "System Patched0");
    ftl_write_full_path(0, host_buffer);
    strcpy((char *)host_buffer, "System Patched1");
    ftl_write_full_path(1, host_buffer);
    strcpy((char *)host_buffer, "System Patched2");
    ftl_write_full_path(2, host_buffer);
    strcpy((char *)host_buffer, "System Patched3");
    ftl_write_full_path(3, host_buffer);
    strcpy((char *)host_buffer, "System Patched4");
    ftl_write_full_path(4, host_buffer);

    // 【新增测试】：直观观察 Vth 分布的劣化过程
    printf("\n--- 可视化阶段：物理 Vth 分布对比 ---\n");

    // 打印刚写完的健康分布 (两座山峰清晰分离)
    print_vth_histogram(g_l2p_table[0], "新盘写入状态 - Health");

    // 暴力老化测试：模拟被疯狂读取 300 次，并伴随时间衰减
    force_age_page(g_l2p_table[0], 200);

    // 打印老化后的分布 (两座山峰变宽，甚至向中心 0.0V 交叠)
    print_vth_histogram(g_l2p_table[0], "严重磨损状态 - Degraded");

    // 3. 【终极毁灭】：Double Program 注入
    force_double_program_page(g_l2p_table[4]);
    print_vth_histogram(g_l2p_table[4], "致命 Bug - Double Program");

    printf("\n--- 第二阶段：模拟磨损后的主机读取验证 ---\n");

    Inject_err_t err = {0};

    // 正常读取刚才被我们暴力老化过的 PPA 0
    ftl_read_full_path(0, err);

    err.inject_erase_error = 1;
    ftl_read_full_path(0, err);
    memset(&err, 0, sizeof(Inject_err_t));
    err.inject_all_zero_error = 1;
    ftl_read_full_path(1, err);
    memset(&err, 0, sizeof(Inject_err_t));
    err.inject_unc_error = 1;
    ftl_read_full_path(2, err);
    memset(&err, 0, sizeof(Inject_err_t));
    err.inject_crc_error = 1;
    ftl_read_full_path(3, err);
    memset(&err, 0, sizeof(Inject_err_t));
    err.inject_meta_mismatch = 1;
    ftl_read_full_path(3, err);

    memset(&err, 0, sizeof(Inject_err_t));
    ftl_read_full_path(4, err);

    printf("\n[Simulator] 全链路物理仿真演示完毕。\n");
}

int main(void)
{
    srand((unsigned)time(NULL));
    printf("========================================================\n");
    printf(" 企业级 SSD 固件全数据通路模拟 (物理 Vth 模型 + LDPC)\n");
    printf("========================================================\n");

    build_mock_rom_table();
    if (ldpc_engine_init() < 0)
        return -1;
    ftl_init();

    ftl_rw_test();

    for (int i = 0; i < PARITY_BITS; i++)
        free(g_check_nodes[i].connected_vns);
    for (int i = 0; i < CODEWORD_BITS; i++)
        free(g_var_nodes[i].connected_cns);
    free(g_check_nodes);
    free(g_var_nodes);
    free(g_rom_base_matrix);

    return 0;
}