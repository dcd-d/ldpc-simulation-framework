#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ==========================================
 * 物理 NAND 与系统级参数定义
 * ========================================== */
// --- 数据通路参数 ---
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

// --- LDPC 矩阵参数 ---
#define LIFTING_FACTOR 64
#define NUM_DATA_BLOCKS (LDPC_PAYLOAD_BITS / LIFTING_FACTOR)
#define NUM_CHK_BLOCKS (PARITY_BITS / LIFTING_FACTOR)
#define COLUMN_WEIGHT 4
#define MAX_ITERATIONS 20
#define MAX_HARD_ITERS 1000

// --- FTL 参数 ---
#define PAGES_PER_BLOCK 64
#define NUM_PHYSICAL_BLOCKS 32
#define NUM_PHYSICAL_PAGES (NUM_PHYSICAL_BLOCKS * PAGES_PER_BLOCK)
#define NUM_LOGICAL_PAGES 1024
#define INVALID_PPA 0xFFFFFFFF

typedef unsigned int LBA;
typedef unsigned int PPA;

/* ==========================================
 * 全局数据结构
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

/* 真实的物理 NAND 阵列 (存放纯物理码字) */
typedef struct
{
    int is_valid;
    unsigned char raw_cell_data[CODEWORD_BYTES];
} PhysicalPage;

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

/* ==========================================
 * 模块 1: 数据完整性引擎 (CRC & Scrambler)
 * ========================================== */
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

void nand_data_randomizer(unsigned char *data, int length, unsigned int seed)
{
    unsigned int lfsr = seed & 0x7FFF;
    lfsr = (seed ^ 0x5AA5) & 0x7FFF; 
    if (lfsr == 0)
        lfsr = 1;
    for (int i = 0; i < length; i++)
    {
        unsigned char random_byte = 0;
        for (int b = 0; b < 8; b++)
        {
            unsigned int bit = ((lfsr >> 14) ^ (lfsr >> 13)) & 1;
            lfsr = ((lfsr << 1) | bit) & 0x7FFF;
            random_byte = (random_byte << 1) | bit;
        }
        data[i] ^= random_byte;
    }
}

/* ==========================================
 * 模块 2: LDPC 引擎核心 (0环矩阵生成 + 编解码)
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

const double g_llr_lut[8] = {15.0, 8.0, 0.0, 3.0, -15.0, -8.0, 0.0, -3.0};
void generate_quantized_llr(const unsigned char *tx_codeword, double *llr_out, double snr_db)
{
    double noise_sigma = 1.0 / sqrt(2.0 * pow(10, snr_db / 10.0));
    for (int i = 0; i < CODEWORD_BITS; i++)
    {
        double bpsk_symbol = (tx_codeword[i] == 0) ? 1.0 : -1.0;
        double awgn_noise = noise_sigma * sqrt(-2.0 * log(((double)rand() / RAND_MAX) + 1e-9)) * cos(2.0 * M_PI * ((double)rand() / RAND_MAX));
        double rx_symbol = bpsk_symbol + awgn_noise;
        unsigned char hb = (rx_symbol < 0.0) ? 1 : 0;
        unsigned char sb1 = (rx_symbol > -0.2 && rx_symbol <= 0.2) ? 1 : 0;
        unsigned char sb2 = (rx_symbol > -0.6 && rx_symbol <= 0.6) ? 1 : 0;
        llr_out[i] = g_llr_lut[(hb << 2) | (sb1 << 1) | sb2];
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
        memset(g_nand_flash[i].raw_cell_data, 0, CODEWORD_BYTES);
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

/* FTL 写入接口：打包 -> 加扰 -> 编码 -> 落盘 */
void ftl_write_full_path(LBA lba, const unsigned char *host_user_data)
{
    PPA new_ppa = allocate_physical_page();
    if (g_l2p_table[lba] != INVALID_PPA)
        g_nand_flash[g_l2p_table[lba]].is_valid = 0; // 标记脏块

    unsigned char ftl_meta_data[META_DATA_BYTES] = {0};
    memcpy(ftl_meta_data, &lba, sizeof(LBA));

    unsigned char payload_bytes[LDPC_PAYLOAD_BYTES];
    pack_nand_page(host_user_data, ftl_meta_data, payload_bytes);

    // 加扰 (使用 PPA 保证唯一性)
    nand_data_randomizer(payload_bytes, LDPC_PAYLOAD_BYTES, new_ppa);

    unsigned char payload_bits[LDPC_PAYLOAD_BITS];
    unsigned char codeword_bits[CODEWORD_BITS];
    bytes_to_bits(payload_bytes, payload_bits, LDPC_PAYLOAD_BYTES);

    ldpc_encode(payload_bits, codeword_bits);
    bits_to_bytes(codeword_bits, g_nand_flash[new_ppa].raw_cell_data, CODEWORD_BYTES);
    g_nand_flash[new_ppa].is_valid = 1;
    g_l2p_table[lba] = new_ppa;

    printf("[FW Write] LBA %04d 映射至 PPA %04d | 净荷封装->LFSR加扰->LDPC编码完毕。\n", lba, new_ppa);
}

/* FTL 读取接口：加噪 -> 解码 -> 解扰 -> 验毒 -> 查映射 */
void ftl_read_full_path(LBA lba)
{
    PPA target_ppa = g_l2p_table[lba];
    if (target_ppa == INVALID_PPA)
    {
        printf("[FW Read]  LBA %04d => [空数据]\n", lba);
        return;
    }

    printf("[FW Read]  请求读取 LBA %04d (物理寻址: PPA %04d)...\n", lba, target_ppa);

    // 1. 模拟信道读取 (设定此时已磨损至 5.0dB)
    double snr_simulate = 5.0;
    unsigned char tx_codeword_bits[CODEWORD_BITS];
    double rx_llr[CODEWORD_BITS];
    unsigned char rx_soft_out_bits[CODEWORD_BITS];

    bytes_to_bits(g_nand_flash[target_ppa].raw_cell_data, tx_codeword_bits, CODEWORD_BYTES);
    generate_quantized_llr(tx_codeword_bits, rx_llr, snr_simulate);

    // 2. LDPC 纠错
    int iter = ldpc_decode_soft(rx_llr, rx_soft_out_bits);
    if (iter < 0)
    {
        printf("  -> [Error] LDPC UE (不可纠错错误)，触发 RAIN 重构！\n");
        return;
    }

    // 3. 提取 Byte 并解扰
    unsigned char decoded_payload[LDPC_PAYLOAD_BYTES];
    bits_to_bytes(rx_soft_out_bits, decoded_payload, LDPC_PAYLOAD_BYTES);
    nand_data_randomizer(decoded_payload, LDPC_PAYLOAD_BYTES, target_ppa);

    // 4. 双重 CRC 验毒
    int crc_res = verify_nand_page(decoded_payload);
    if (crc_res == -1)
    {
        printf("  -> [Error] User CRC 校验失败 (发生了静默数据损坏 SDC)！\n");
        return;
    }
    if (crc_res == -2)
    {
        printf("  -> [Error] Meta CRC 校验失败！\n");
        return;
    }

    // 5. 校验 FTL 元数据反向映射
    LBA read_lba;
    memcpy(&read_lba, decoded_payload + USER_DATA_BYTES + USER_CRC_BYTES, sizeof(LBA));
    if (read_lba != lba)
    {
        printf("  -> [Error] 地址错乱！读出的 OOB 指向 LBA %d！\n", read_lba);
        return;
    }

    // 6. 成功提纯 User Data
    unsigned char final_host_data[16]; // 这里只截取前15个字符打印演示
    memcpy(final_host_data, decoded_payload, 15);
    final_host_data[15] = '\0';

    printf("  -> [Success] 纠错耗时 %d 轮。CRC 验毒通过！提取前缀数据: \"%s\"\n", iter, final_host_data);
}

/* ==========================================
 * 测试引擎：体验完整的 SSD 固件工作流
 * ========================================== */
int main(void)
{
    srand((unsigned)time(NULL));

    printf("========================================================\n");
    printf(" 企业级 SSD 固件全数据通路模拟 (FTL + LDPC + PHY)\n");
    printf("========================================================\n");

    /* 初始化底层图模型与 FTL */
    build_mock_rom_table();
    if (ldpc_engine_init() < 0)
        return -1;
    ftl_init();

    printf("\n--- 第一阶段：主机初始化写入 ---\n");
    unsigned char host_buffer[USER_DATA_BYTES] = {0};

    strcpy((char *)host_buffer, "System Boot Up!");
    ftl_write_full_path(10, host_buffer);

    strcpy((char *)host_buffer, "File Header 100");
    ftl_write_full_path(25, host_buffer);

    printf("\n--- 第二阶段：异地更新 (Out-of-Place Update) ---\n");
    strcpy((char *)host_buffer, "System Patched.");
    ftl_write_full_path(10, host_buffer);

    printf("\n--- 第三阶段：模拟磨损后的主机读取验证 ---\n");
    ftl_read_full_path(25);
    ftl_read_full_path(10);
    ftl_read_full_path(99); // 读空盘

    printf("\n[Simulator] 全链路演示完毕。\n");

    /* 清理内存 */
    for (int i = 0; i < PARITY_BITS; i++)
        free(g_check_nodes[i].connected_vns);
    for (int i = 0; i < CODEWORD_BITS; i++)
        free(g_var_nodes[i].connected_cns);
    free(g_check_nodes);
    free(g_var_nodes);
    free(g_rom_base_matrix);

    return 0;
}