#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ==========================================
 * 物理参数宏定义 (推荐 Debug 配置)
 * ========================================== */
// #define DATA_BYTES 8
#define DATA_BYTES 4096

#if DATA_BYTES == 4096
#define DATA_BYTES 4096
#define PARITY_BYTES 512
#define DATA_BITS (DATA_BYTES * 8)
#define PARITY_BITS (PARITY_BYTES * 8)
#define CODEWORD_BITS (DATA_BITS + PARITY_BITS)

#define LIFTING_FACTOR 64
#define NUM_DATA_BLOCKS (DATA_BITS / LIFTING_FACTOR)
#define NUM_CHK_BLOCKS (PARITY_BITS / LIFTING_FACTOR)
#define COLUMN_WEIGHT 4

#define MAX_ITERATIONS 20
#define MAX_HARD_ITERS 1000

#elif DATA_BYTES == 8
/* ==========================================
 * 推荐 Debug 配置 (Data: 8B, Parity: 2B)
 * ========================================== */
#define DATA_BYTES 8
#define PARITY_BYTES 2

#define DATA_BITS (DATA_BYTES * 8)              // 64 bits
#define PARITY_BITS (PARITY_BYTES * 8)          // 16 bits
#define CODEWORD_BITS (DATA_BITS + PARITY_BITS) // 80 bits

#define LIFTING_FACTOR 4                              // 矩阵由 4x4 的小块组成
#define NUM_DATA_BLOCKS (DATA_BITS / LIFTING_FACTOR)  // 16 列
#define NUM_CHK_BLOCKS (PARITY_BITS / LIFTING_FACTOR) // 4 行
#define COLUMN_WEIGHT 3                               // 每列 3 个非零块

#define MAX_ITERATIONS 100
#define MAX_HARD_ITERS 5000
#endif

/* ==========================================
 * 数据结构与全局状态
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

/* 硬件引擎 SRAM 状态与 ROM 表 */
static struct ldpc_rom_entry *g_rom_base_matrix = NULL;
static int g_rom_matrix_entries = 0;
static struct check_node *g_check_nodes = NULL;
static struct var_node *g_var_nodes = NULL;
static int g_is_initialized = 0;

/* ==========================================
 * 模块 1: 环境初始化与内存管理
 * ========================================== */
void build_mock_rom_table(void)
{
    int j, w, row_idx;
    int total_entries = NUM_DATA_BLOCKS * COLUMN_WEIGHT;
    int current_entry = 0;
    int temp_matrix[NUM_CHK_BLOCKS][NUM_DATA_BLOCKS];

    g_rom_base_matrix = (struct ldpc_rom_entry *)malloc(
        total_entries * sizeof(struct ldpc_rom_entry));
    memset(temp_matrix, -1, sizeof(temp_matrix));

    // srand(12345); // 固定随机种子，保证每次调试矩阵结构一致
    for (j = 0; j < NUM_DATA_BLOCKS; j++)
    {
        for (w = 0; w < COLUMN_WEIGHT; w++)
        {
            row_idx = rand() % NUM_CHK_BLOCKS;
            if (temp_matrix[row_idx][j] == -1)
            {
                temp_matrix[row_idx][j] = rand() % LIFTING_FACTOR;
                g_rom_base_matrix[current_entry].row = row_idx;
                g_rom_base_matrix[current_entry].col = j;
                g_rom_base_matrix[current_entry].shift = temp_matrix[row_idx][j];
                current_entry++;
            }
            else
            {
                w--; /* 重新找一个空行 */
            }
        }
    }
    g_rom_matrix_entries = current_entry;
    // printf("[Testbed] 成功生成模拟 ROM 表，包含 %d 个有效宏块。\n",
    //        g_rom_matrix_entries);
}

/* ==========================================
 * 模块 1: 工业级 QC-LDPC 无 4 环矩阵生成器
 * ========================================== */
void build_mock_rom_table_qc(void)
{
    int j, w, row_idx, shift_val;
    int c_prime, r_prime, diff;
    int total_entries = NUM_DATA_BLOCKS * COLUMN_WEIGHT;
    int current_entry = 0;
    int temp_matrix[NUM_CHK_BLOCKS][NUM_DATA_BLOCKS];
    int is_valid, attempts, restart_count = 0;

restart_matrix:
    if (restart_count > 5000) {
        printf("[Fatal] Z = %d 空间太小，无法在当前列重下完全避开 4 环！\n", LIFTING_FACTOR);
        exit(-1);
    }

    memset(temp_matrix, -1, sizeof(temp_matrix));
    current_entry = 0;

    for (j = 0; j < NUM_DATA_BLOCKS; j++)
    {
        for (w = 0; w < COLUMN_WEIGHT; w++)
        {
            /* 1. 随机找一个还没填过数据的空行 */
            int row_attempts = 0;
            do {
                row_idx = rand() % NUM_CHK_BLOCKS;
                row_attempts++;
                // 如果随机 100 次都找不到空位，说明这列塞满了，死锁重启
                if (row_attempts > 100) goto restart_matrix; 
            } while (temp_matrix[row_idx][j] != -1);

            /* 2. 尝试分配 Shift 移位值，并进行严格的代数同余校验 */
            attempts = 0;
            do {
                shift_val = rand() % LIFTING_FACTOR;
                is_valid = 1;

                // 遍历之前所有的列 (c_prime)
                for (c_prime = 0; c_prime < j; c_prime++) {
                    // 必须在这两列的同一行都有数据，才可能构成 2x2 网格
                    if (temp_matrix[row_idx][c_prime] != -1) {
                        // 寻找另一个同时连接这两列的行 (r_prime)
                        for (r_prime = 0; r_prime < NUM_CHK_BLOCKS; r_prime++) {
                            if (r_prime != row_idx &&
                                temp_matrix[r_prime][j] != -1 &&
                                temp_matrix[r_prime][c_prime] != -1)
                            {
                                // 抓到了一个 2x2 的宏块！提取它们的 Shift 值
                                int s11 = temp_matrix[r_prime][c_prime];
                                int s12 = temp_matrix[r_prime][j];
                                int s21 = temp_matrix[row_idx][c_prime];
                                int s22 = shift_val; 

                                // 核心公式: (S11 + S22) - (S12 + S21) = 0 mod Z
                                diff = (s11 + s22) - (s12 + s21);
                                if ((diff % LIFTING_FACTOR + LIFTING_FACTOR) % LIFTING_FACTOR == 0) {
                                    is_valid = 0; // 产生 4 环！一票否决
                                    break;
                                }
                            }
                        }
                    }
                    if (!is_valid) break;
                }
                
                attempts++;
                // 如果试遍了所有的 Shift 值都不行，说明走进了死胡同，全盘掀桌子重来！
                if (attempts > LIFTING_FACTOR * 2) {
                    restart_count++;
                    goto restart_matrix;
                }
            } while (!is_valid);

            /* 3. 代数校验通过，落子无悔 */
            temp_matrix[row_idx][j] = shift_val;
        }
    }

    /* 将二维数组转换为一维 ROM 表格给硬件引擎使用 */
    g_rom_base_matrix = (struct ldpc_rom_entry *)malloc(total_entries * sizeof(struct ldpc_rom_entry));
    for (j = 0; j < NUM_DATA_BLOCKS; j++) {
        for (int i = 0; i < NUM_CHK_BLOCKS; i++) {
            if (temp_matrix[i][j] != -1) {
                g_rom_base_matrix[current_entry].row = i;
                g_rom_base_matrix[current_entry].col = j;
                g_rom_base_matrix[current_entry].shift = temp_matrix[i][j];
                current_entry++;
            }
        }
    }
    g_rom_matrix_entries = current_entry;
    printf("[Matrix Generator] 历经 %d 次死锁回溯，成功生成【0 个 4 环】的完美 QC 矩阵！\n", restart_count);
}

int ldpc_engine_init(void)
{
    int i, entry_idx, cn_idx, vn_idx, z_offset;

    if (g_is_initialized)
        return 0;

    g_check_nodes =
        (struct check_node *)calloc(PARITY_BITS, sizeof(struct check_node));
    g_var_nodes =
        (struct var_node *)calloc(CODEWORD_BITS, sizeof(struct var_node));

    if (!g_check_nodes || !g_var_nodes)
        return -1;

    /* 第一遍扫描：统计度数 */
    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            cn_idx = g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR + z_offset;
            vn_idx = g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR +
                     (z_offset + g_rom_base_matrix[entry_idx].shift) % LIFTING_FACTOR;
            g_check_nodes[cn_idx].degree++;
            g_var_nodes[vn_idx].degree++;
        }
    }

    /* 校验区对角线度数补充 */
    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].degree++;
        g_var_nodes[DATA_BITS + i].degree++;
    }

    /* 精准分配内存并重置游标 */
    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].connected_vns =
            (int *)malloc(g_check_nodes[i].degree * sizeof(int));
        g_check_nodes[i].degree = 0;
    }
    for (i = 0; i < CODEWORD_BITS; i++)
    {
        g_var_nodes[i].connected_cns =
            (int *)malloc(g_var_nodes[i].degree * sizeof(int));
        g_var_nodes[i].degree = 0;
    }

    /* 第二遍扫描：双向连线 */
    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            cn_idx = g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR + z_offset;
            vn_idx = g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR +
                     (z_offset + g_rom_base_matrix[entry_idx].shift) % LIFTING_FACTOR;

            g_check_nodes[cn_idx].connected_vns[g_check_nodes[cn_idx].degree++] =
                vn_idx;
            g_var_nodes[vn_idx].connected_cns[g_var_nodes[vn_idx].degree++] = cn_idx;
        }
    }

    /* 校验区对角线连线补充 */
    for (i = 0; i < PARITY_BITS; i++)
    {
        g_check_nodes[i].connected_vns[g_check_nodes[i].degree++] = DATA_BITS + i;
        g_var_nodes[DATA_BITS + i]
            .connected_cns[g_var_nodes[DATA_BITS + i].degree++] = i;
    }

    g_is_initialized = 1;
    // printf("[LDPC Engine] 硬件图模型初始化完毕。\n");
    return 0;
}

void ldpc_engine_cleanup(void)
{
    int i;
    if (!g_is_initialized)
        return;

    for (i = 0; i < PARITY_BITS; i++)
        free(g_check_nodes[i].connected_vns);
    for (i = 0; i < CODEWORD_BITS; i++)
        free(g_var_nodes[i].connected_cns);
    free(g_check_nodes);
    free(g_var_nodes);
    if (g_rom_base_matrix)
        free(g_rom_base_matrix);
    g_is_initialized = 0;
}

/* ==========================================
 * 模块 2: 编码与信道模拟
 * ========================================== */
void ldpc_encode(const unsigned char *data_in, unsigned char *codeword_out)
{
    int entry_idx, z_offset, data_bit_idx, parity_bit_idx;

    memcpy(codeword_out, data_in, DATA_BITS);
    memset(codeword_out + DATA_BITS, 0, PARITY_BITS);

    for (entry_idx = 0; entry_idx < g_rom_matrix_entries; entry_idx++)
    {
        for (z_offset = 0; z_offset < LIFTING_FACTOR; z_offset++)
        {
            data_bit_idx =
                g_rom_base_matrix[entry_idx].col * LIFTING_FACTOR + z_offset;
            parity_bit_idx =
                g_rom_base_matrix[entry_idx].row * LIFTING_FACTOR +
                (z_offset + LIFTING_FACTOR - g_rom_base_matrix[entry_idx].shift) %
                    LIFTING_FACTOR;

            codeword_out[DATA_BITS + parity_bit_idx] ^= data_in[data_bit_idx];
        }
    }
}

void generate_llr(const unsigned char *tx_codeword, double *llr_out,
                  double snr_db)
{
    int i;
    double snr_linear = pow(10, snr_db / 10.0);
    double noise_sigma = 1.0 / sqrt(2.0 * snr_linear);
    double bpsk_symbol, rand_u1, rand_u2, awgn_noise, rx_symbol;

    for (i = 0; i < CODEWORD_BITS; i++)
    {
        bpsk_symbol = (tx_codeword[i] == 0) ? 1.0 : -1.0;
        /* Box-Muller 变换生成高斯白噪声 */
        rand_u1 = (double)rand() / RAND_MAX;
        rand_u2 = (double)rand() / RAND_MAX;
        awgn_noise = noise_sigma * sqrt(-2.0 * log(rand_u1 + 1e-9)) *
                     cos(2.0 * M_PI * rand_u2);

        rx_symbol = bpsk_symbol + awgn_noise;
        llr_out[i] = 2.0 * rx_symbol / (noise_sigma * noise_sigma);
    }
}

/* ==========================================
 * 真实的 NAND 3-bit 软读取查表 (Look-Up Table)
 * ========================================== */
// 索引: {HB, SB1, SB2} 的 3-bit 组合
const double g_llr_lut[8] = {
    15.0,   // 000: 极强 0 (电压 > 0.6)
    8.0,    // 001: 中等 0 (0.2 < 电压 <= 0.6)
    0.0,    // 010: 非法状态 (物理上不可能出现的极性)
    3.0,    // 011: 弱 0   (0 < 电压 <= 0.2)
    -15.0,  // 100: 极强 1 (电压 < -0.6)
    -8.0,   // 101: 中等 1 (-0.6 <= 电压 < -0.2)
    0.0,    // 110: 非法状态
    -3.0    // 111: 弱 1   (-0.2 <= 电压 <= 0)
};

void generate_quantized_llr(const unsigned char *tx_codeword, double *llr_out, double snr_db)
{
    int i;
    double snr_linear = pow(10, snr_db / 10.0);
    double noise_sigma = 1.0 / sqrt(2.0 * snr_linear);
    double bpsk_symbol, rand_u1, rand_u2, awgn_noise, rx_symbol;
    unsigned char hb, sb1, sb2, lut_index;

    for (i = 0; i < CODEWORD_BITS; i++)
    {
        /* 1. 物理层：模拟真实的细胞电压 (带噪声) */
        bpsk_symbol = (tx_codeword[i] == 0) ? 1.0 : -1.0;
        rand_u1 = (double)rand() / RAND_MAX;
        rand_u2 = (double)rand() / RAND_MAX;
        awgn_noise = noise_sigma * sqrt(-2.0 * log(rand_u1 + 1e-9)) * cos(2.0 * M_PI * rand_u2);
        rx_symbol = bpsk_symbol + awgn_noise;

        /* 2. NAND 接口层：模拟 3-Step Soft Read (电压切片) */
        // 假设中心读取电压 V0 = 0，内圈读取 V_inner = ±0.2，外圈读取 V_outer = ±0.6
        
        // HB (Hard Bit): 以 0 为界限
        hb = (rx_symbol < 0.0) ? 1 : 0;
        
        // SB1 (Soft Bit 1 - 内圈): 落在 [-0.2, 0.2] 的最模糊地带为 1，否则为 0
        sb1 = (rx_symbol > -0.2 && rx_symbol <= 0.2) ? 1 : 0;
        
        // SB2 (Soft Bit 2 - 外圈): 落在 [-0.6, 0.6] 之间为 1，否则为 0
        sb2 = (rx_symbol > -0.6 && rx_symbol <= 0.6) ? 1 : 0;

        /* 3. 固件逻辑层：位运算拼凑 3-bit 索引 */
        lut_index = (hb << 2) | (sb1 << 1) | sb2;

        /* 4. 译码器硬件层：查表输出离散 LLR！彻底抛弃公式！ */
        llr_out[i] = g_llr_lut[lut_index];
    }
}

int ldpc_decode_hard(unsigned char *rx_codeword)
{
    #if 0
    ldpc_decode_hard_pbf(rx_codeword);
    #else
    ldpc_decode_hard_old(rx_codeword);
    #endif
}
/* ==========================================
 * 模块 3: 核心解码算法 (升级版 PBF 硬判决)
 * ========================================== */
int ldpc_decode_hard_pbf(unsigned char *rx_codeword)
{
    unsigned char *syndrome_buf = (unsigned char *)malloc(PARITY_BITS);
    int *fail_counts = (int *)calloc(CODEWORD_BITS, sizeof(int));
    int iter, i, edge_idx, is_all_zero, max_fail_peak;
    unsigned char current_syndrome;
    
    // --- 新增：死锁追踪器 ---
    int last_syndrome_weight = 1e9; // 记录上一轮有多少个方程报错
    int deadlock_counter = 0;       // 记录原地踏步了多少轮

    for (iter = 0; iter < MAX_HARD_ITERS; iter++)
    {
        is_all_zero = 1;
        int current_syndrome_weight = 0; // 统计当前轮的报错总数
        
        /* 1. 计算伴随式 */
        for (i = 0; i < PARITY_BITS; i++)
        {
            current_syndrome = 0;
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                current_syndrome ^= rx_codeword[g_check_nodes[i].connected_vns[edge_idx]];
            }
            syndrome_buf[i] = current_syndrome;
            if (current_syndrome != 0) {
                is_all_zero = 0;
                current_syndrome_weight++; // 记录报错总数
            }
        }

        if (is_all_zero)
        {
            free(syndrome_buf);
            free(fail_counts);
            return iter + 1;
        }

        /* --- 新增：死锁状态评估 --- */
        // 如果当前报错的方程数量，没有比上一轮少，说明我们卡在陷阱集里了
        if (current_syndrome_weight >= last_syndrome_weight) {
            deadlock_counter++;
        } else {
            deadlock_counter = 0; // 一旦有进展，立刻清零死锁计数
        }
        last_syndrome_weight = current_syndrome_weight;

        /* 2. 统计变量节点报错次数 */
        max_fail_peak = 0;
        memset(fail_counts, 0, CODEWORD_BITS * sizeof(int));

        for (i = 0; i < CODEWORD_BITS; i++)
        {
            for (edge_idx = 0; edge_idx < g_var_nodes[i].degree; edge_idx++)
            {
                if (syndrome_buf[g_var_nodes[i].connected_cns[edge_idx]])
                    fail_counts[i]++;
            }
            if (fail_counts[i] > max_fail_peak)
                max_fail_peak = fail_counts[i];
        }

        /* 3. 智能比特翻转 (引入概率打破对称性) */
        if (max_fail_peak > 0)
        {
            for (i = 0; i < CODEWORD_BITS; i++)
            {
                if (fail_counts[i] == max_fail_peak)
                {
                    // 如果发现连续震荡超过 2 轮，我们就不再 100% 翻转
                    // 而是用 80% 的概率去翻转它，强行撕裂陷阱集的对称结构
                    if (deadlock_counter > 2) {
                        if (rand() % 100 < 80) {
                            rx_codeword[i] ^= 1;
                        }
                    } else {
                        // 正常情况下，保持贪心的高效翻转
                        rx_codeword[i] ^= 1;
                    }
                }
            }
        }
    }

    free(syndrome_buf);
    free(fail_counts);
    return -1;
}
/* ==========================================
 * 模块 3: 核心解码算法
 * ========================================== */
int ldpc_decode_hard_old(unsigned char *rx_codeword)
{
    unsigned char *syndrome_buf = (unsigned char *)malloc(PARITY_BITS);
    int *fail_counts = (int *)calloc(CODEWORD_BITS, sizeof(int));
    int iter, i, edge_idx, is_all_zero, max_fail_peak;
    unsigned char current_syndrome;

    for (iter = 0; iter < MAX_HARD_ITERS; iter++)
    {
        is_all_zero = 1;
        /* 1. 计算伴随式 */
        for (i = 0; i < PARITY_BITS; i++)
        {
            current_syndrome = 0;
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                current_syndrome ^=
                    rx_codeword[g_check_nodes[i].connected_vns[edge_idx]];
            }
            syndrome_buf[i] = current_syndrome;
            if (current_syndrome != 0)
                is_all_zero = 0;
        }

        if (is_all_zero)
        {
            free(syndrome_buf);
            free(fail_counts);
            return iter + 1;
        }

        /* 2. 统计变量节点报错次数 */
        max_fail_peak = 0;
        memset(fail_counts, 0, CODEWORD_BITS * sizeof(int));

        for (i = 0; i < CODEWORD_BITS; i++)
        {
            for (edge_idx = 0; edge_idx < g_var_nodes[i].degree; edge_idx++)
            {
                if (syndrome_buf[g_var_nodes[i].connected_cns[edge_idx]])
                    fail_counts[i]++;
            }
            if (fail_counts[i] > max_fail_peak)
                max_fail_peak = fail_counts[i];
        }

        /* 3. 翻转错误比特 */
        if (max_fail_peak > 0)
        {
            for (i = 0; i < CODEWORD_BITS; i++)
            {
                if (fail_counts[i] == max_fail_peak)
                    rx_codeword[i] ^= 1;
            }
        }
    }

    free(syndrome_buf);
    free(fail_counts);
    return -1;
}

static double extract_sign(double value) { return (value >= 0) ? 1.0 : -1.0; }

int ldpc_decode_soft(double *channel_llr, unsigned char *out_codeword)
{
    double **msg_chk_to_var;
    double *posteriori_llr, *msg_var_to_chk;
    int iter, i, j, edge_idx, search_idx;
    int max_node_degree = 0, is_all_zero;
    int target_vn, target_cn, back_edge_idx, min_mag1_idx;
    double min_mag1, min_mag2, total_sign_prod;
    double current_mag, excl_sign, excl_mag;
    unsigned char current_syndrome;

    msg_chk_to_var = (double **)malloc(PARITY_BITS * sizeof(double *));
    for (i = 0; i < PARITY_BITS; i++)
    {
        msg_chk_to_var[i] =
            (double *)calloc(g_check_nodes[i].degree, sizeof(double));
        if (g_check_nodes[i].degree > max_node_degree)
            max_node_degree = g_check_nodes[i].degree;
    }

    posteriori_llr = (double *)malloc(CODEWORD_BITS * sizeof(double));
    msg_var_to_chk = (double *)malloc(max_node_degree * sizeof(double));

    for (iter = 0; iter < MAX_ITERATIONS; iter++)
    {
        /* 1. 校验节点更新 (Row Processing) */
        for (i = 0; i < PARITY_BITS; i++)
        {
            min_mag1 = 1e9;
            min_mag2 = 1e9;
            min_mag1_idx = -1;
            total_sign_prod = 1.0;

            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                target_vn = g_check_nodes[i].connected_vns[edge_idx];
                msg_var_to_chk[edge_idx] =
                    (iter == 0)
                        ? channel_llr[target_vn]
                        : (posteriori_llr[target_vn] - msg_chk_to_var[i][edge_idx]);

                current_mag = fabs(msg_var_to_chk[edge_idx]);
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
                excl_sign = total_sign_prod * extract_sign(msg_var_to_chk[edge_idx]);
                excl_mag = (edge_idx == min_mag1_idx) ? min_mag2 : min_mag1;
                /* 归一化 Min-Sum 缩放因子 0.75 */
                msg_chk_to_var[i][edge_idx] = 0.75 * excl_sign * excl_mag;
            }
        }

        /* 2. 变量节点更新 (Column Processing) */
        is_all_zero = 1;
        for (j = 0; j < CODEWORD_BITS; j++)
        {
            posteriori_llr[j] = channel_llr[j];
            for (edge_idx = 0; edge_idx < g_var_nodes[j].degree; edge_idx++)
            {
                target_cn = g_var_nodes[j].connected_cns[edge_idx];
                back_edge_idx = -1;
                for (search_idx = 0; search_idx < g_check_nodes[target_cn].degree;
                     search_idx++)
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

        /* 3. 伴随式检查提前退出 */
        for (i = 0; i < PARITY_BITS; i++)
        {
            current_syndrome = 0;
            for (edge_idx = 0; edge_idx < g_check_nodes[i].degree; edge_idx++)
            {
                current_syndrome ^=
                    out_codeword[g_check_nodes[i].connected_vns[edge_idx]];
            }
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
 * 模块 3.5: 矩阵拓扑优化器 (消除 Cycle-4)
 * ========================================== */

// 扫描当前图模型中的 4 环数量
int count_cycle_4(void)
{
    int cycle4_count = 0;
    int i, e1, e2, c1, c2, r, c;
    
    /* 分配一维连续内存用于统计 (4096*4096*2 bytes ≈ 33.5MB) */
    unsigned short *shared_vns = (unsigned short *)calloc(PARITY_BITS * PARITY_BITS, sizeof(unsigned short));

    /* 遍历所有变量节点，两两组合它连接的校验节点 */
    for (i = 0; i < CODEWORD_BITS; i++) {
        for (e1 = 0; e1 < g_var_nodes[i].degree; e1++) {
            c1 = g_var_nodes[i].connected_cns[e1];
            for (e2 = e1 + 1; e2 < g_var_nodes[i].degree; e2++) {
                c2 = g_var_nodes[i].connected_cns[e2];
                // 确保行坐标小于列坐标，只用矩阵的上三角
                r = (c1 < c2) ? c1 : c2;
                c = (c1 > c2) ? c1 : c2;
                // 记录这两个校验节点共同连接的变量节点数
                shared_vns[r * PARITY_BITS + c]++;
            }
        }
    }

    /* 如果有两个以上的变量节点连向了相同的一对校验节点，就算作 4 环 */
    for (i = 0; i < PARITY_BITS * PARITY_BITS; i++) {
        if (shared_vns[i] > 1) {
            int k = shared_vns[i];
            // 从 k 个节点中选 2 个的组合数: C(k, 2)
            cycle4_count += (k * (k - 1)) / 2; 
        }
    }

    free(shared_vns);
    return cycle4_count;
}

// 蒙特卡洛搜索最优矩阵
void find_optimal_matrix(int search_trials)
{
    int best_cycle4 = 2e9; // 初始设为无限大
    struct ldpc_rom_entry *best_matrix = NULL;
    int best_entries_count = 0;
    int t, current_cycle4;

    printf("\n[Matrix Optimizer] 开始搜索最优 ROM 矩阵 (基于 Cycle-4 最小化)，尝试 %d 次...\n", search_trials);

    for (t = 0; t < search_trials; t++)
    {
        if (t % 100 == 0)
        {
            printf("serch %d times\n", t);
        }
        // 1. 生成一版随机矩阵
        build_mock_rom_table();
        
        // 2. 根据这版矩阵构建物理连线图
        if (ldpc_engine_init() < 0) return;

        // 3. 统计 4 环的致命数量
        current_cycle4 = count_cycle_4();

        // 4. 如果发现更优秀的拓扑结构，保存它！
        if (current_cycle4 < best_cycle4)
        {
            best_cycle4 = current_cycle4;
            if (best_matrix) free(best_matrix);
            
            best_matrix = (struct ldpc_rom_entry *)malloc(g_rom_matrix_entries * sizeof(struct ldpc_rom_entry));
            memcpy(best_matrix, g_rom_base_matrix, g_rom_matrix_entries * sizeof(struct ldpc_rom_entry));
            best_entries_count = g_rom_matrix_entries;
            
            printf("  -> 发现优质矩阵! 迭代 %3d, Cycle-4 数量降低至: %d\n", t, best_cycle4);
        }

        // 清理当前连线，为下一次尝试腾出空间
        ldpc_engine_cleanup();
    }

    printf("[Matrix Optimizer] 搜索完成！选定包含 %d 个 Cycle-4 的矩阵作为固件基线。\n", best_cycle4);

    // 5. 将万里挑一的最优矩阵恢复到全局状态，并正式初始化
    g_rom_base_matrix = best_matrix;
    g_rom_matrix_entries = best_entries_count;
    ldpc_engine_init();
}

/* ==========================================
 * 模块 4: 纠错极限探测 (Stress Test)
 * ========================================== */
void test_correction_limit(const unsigned char *tx_codeword)
{
    unsigned char *hard_rx = (unsigned char *)malloc(CODEWORD_BITS);
    double *rx_llr_channel = (double *)malloc(CODEWORD_BITS * sizeof(double));
    unsigned char *soft_out = (unsigned char *)malloc(CODEWORD_BITS);
    
    int t, i;
    const int trials_per_snr = 20; // 每个 SNR 测试 20 次
    double snr;

    printf("\n[Stress Test] 启动 SNR 扫频测试 (模拟 NAND 磨损加深)...\n");
    printf(" SNR(dB) | 硬判决成功率 | 软判决成功率 | 平均原始误码(Bits)\n");
    printf("----------------------------------------------------------\n");

    // 从极其健康的 6.5dB，一路折磨到极度衰老的 3.5dB
    for (snr = 8.0; snr >= 1.0; snr -= 0.5)
    {
        int hard_ok = 0, soft_ok = 0;
        int total_raw_errs = 0;

        for (t = 0; t < trials_per_snr; t++)
        {
            /* 1. 统一生成带噪声的信道模拟信息 (Soft LLR) */
            // generate_llr(tx_codeword, rx_llr_channel, snr);
            generate_quantized_llr(tx_codeword, rx_llr_channel, snr);

            /* 2. 模拟 Hard Read：基于阈值 0 强行切分出 0 和 1 */
            int raw_errs = 0;
            for (i = 0; i < CODEWORD_BITS; i++) {
                hard_rx[i] = (rx_llr_channel[i] < 0) ? 1 : 0;
                // 顺便统计一下，这个 SNR 下物理层到底碎了多少个 bit
                if (hard_rx[i] != tx_codeword[i]) raw_errs++; 
            }
            total_raw_errs += raw_errs;

            /* 3. 公平竞技开始 */
            // 硬判决只拿到粗糙的 hard_rx 盲解
            if (ldpc_decode_hard(hard_rx) >= 0) hard_ok++;
            
            // 软判决拿到包含丰富犹豫度信息的 LLR 细解
            if (ldpc_decode_soft(rx_llr_channel, soft_out) >= 0) soft_ok++;
        }

        printf("  %4.1f   |    %6.1f%%    |    %6.1f%%    |   ~ %d bits\n", 
               snr, 
               (double)hard_ok/trials_per_snr*100, 
               (double)soft_ok/trials_per_snr*100,
               total_raw_errs / trials_per_snr);

        // 如果两个都彻底死透了，停止扫频
        if (hard_ok == 0 && soft_ok == 0) break;
    }

    printf("----------------------------------------------------------\n");
    free(hard_rx);
    free(rx_llr_channel);
    free(soft_out);
}

/* ==========================================
 * 模块 5: 固件级跨 Die RAID 容灾恢复 (RAIN)
 * ========================================== */
#define NUM_DIES 5  // 4 个数据 Die + 1 个校验 Die

void firmware_raid_recovery_demo(void)
{
    unsigned char *raw_data[NUM_DIES];
    unsigned char *tx_codeword[NUM_DIES];
    double *rx_llr[NUM_DIES];
    unsigned char *rx_soft_out[NUM_DIES];
    unsigned char *recovered_die0 = (unsigned char *)malloc(CODEWORD_BITS);
    int i, d, iter;
    int raid_success = 1;

    printf("\n========================================================\n");
    printf(" [FW RAID Flow] 启动跨 Die 异或容灾恢复演练 (4 Data + 1 Parity)\n");
    printf("========================================================\n");

    /* 1. 初始化 5 个 Die 的内存，并生成用户数据 */
    for (d = 0; d < NUM_DIES; d++) {
        raw_data[d] = (unsigned char *)malloc(DATA_BITS);
        tx_codeword[d] = (unsigned char *)malloc(CODEWORD_BITS);
        rx_llr[d] = (double *)malloc(CODEWORD_BITS * sizeof(double));
        rx_soft_out[d] = (unsigned char *)malloc(CODEWORD_BITS);
    }

    // 填充 Die 0~3 的用户数据并进行 LDPC 编码
    for (d = 0; d < NUM_DIES - 1; d++) {
        for (i = 0; i < DATA_BITS; i++) raw_data[d][i] = rand() % 2;
        ldpc_encode(raw_data[d], tx_codeword[d]);
    }

    // 计算 Die 4 的 RAID 校验页 (按位异或前 4 个 Die 的码字)
    for (i = 0; i < CODEWORD_BITS; i++) {
        tx_codeword[4][i] = tx_codeword[0][i] ^ tx_codeword[1][i] ^ 
                            tx_codeword[2][i] ^ tx_codeword[3][i];
    }
    printf("[RAID Engine] 5 个 Die 数据已条带化写入完毕。\n");

    /* 2. 模拟突发物理灾难 */
    printf("\n[Channel] 模拟岁月侵蚀... Die 0 遭遇极度恶劣的漏电损坏！\n");
    
    // Die 0 给予致死量的噪声 (2.5dB)，绝对超出 LDPC 极限
    generate_quantized_llr(tx_codeword[0], rx_llr[0], 2.5); 
    
    // Die 1~4 给予正常的晚期噪声 (5.5dB)，在 LDPC 软解能力范围内
    for (d = 1; d < NUM_DIES; d++) {
        generate_quantized_llr(tx_codeword[d], rx_llr[d], 5.5);
    }

    /* 3. Host 请求读取 Die 0 */
    printf("\n[Host Req] 主机发起读取 LBA (映射至 Die 0)...\n");
    iter = ldpc_decode_soft(rx_llr[0], rx_soft_out[0]);
    if (iter < 0) {
        printf("[LDPC Core] 警告：Die 0 软判决解码彻底失败 (UE)！\n");
    } else {
        printf("[LDPC Core] Die 0 解码成功 (不可能发生)。\n");
        return;
    }

    /* 4. 固件触发 RAID 恢复状态机 */
    printf("\n[FW Exception] 挂起 Host IO，触发后台 RAID 数据恢复流程...\n");
    for (d = 1; d < NUM_DIES; d++) {
        printf("  -> 正在并行抢读 Die %d ... ", d);
        iter = ldpc_decode_soft(rx_llr[d], rx_soft_out[d]);
        if (iter >= 0) {
            printf("LDPC 软解成功 (耗费 %d 轮)\n", iter);
        } else {
            printf("LDPC 软解失败！\n");
            raid_success = 0;
            break;
        }
    }

    /* 5. 执行异或重构 */
    if (raid_success) {
        printf("\n[RAID Engine] 所有辅助 Die 读取成功，开始执行 XOR 异或重构...\n");
        for (i = 0; i < CODEWORD_BITS; i++) {
            // Die0 = Die1 ^ Die2 ^ Die3 ^ Parity(Die4)
            recovered_die0[i] = rx_soft_out[1][i] ^ rx_soft_out[2][i] ^ 
                                rx_soft_out[3][i] ^ rx_soft_out[4][i];
        }

        // 校验恢复出的数据是否和最初存入的一模一样
        int is_perfect = 1;
        for (i = 0; i < DATA_BITS; i++) {
            if (recovered_die0[i] != raw_data[0][i]) {
                is_perfect = 0; break;
            }
        }

        if (is_perfect) {
            printf("[FW Output] 【神迹重现】Die 0 数据完美恢复！安全返回给 Host！避免了蓝屏灾难。\n");
        } else {
            printf("[FW Output] 恢复失败，数据不一致。\n");
        }
    } else {
        printf("[FW Output] 【彻底死机】超过两个 Die 损坏，RAID 5 阵列崩溃，上报蓝屏 (BSOD)。\n");
    }

    /* 释放内存 */
    for (d = 0; d < NUM_DIES; d++) {
        free(raw_data[d]); free(tx_codeword[d]); 
        free(rx_llr[d]); free(rx_soft_out[d]);
    }
    free(recovered_die0);
}

/* ==========================================
 * 主函数入口
 * ========================================== */
int main(void)
{
    unsigned char *raw_data = (unsigned char *)malloc(DATA_BITS);
    unsigned char *tx_codeword = (unsigned char *)malloc(CODEWORD_BITS);
    int i;

    /* 必须把随机种子放到最开头，让优化器每次都能出不同结果 */
    srand((unsigned)time(NULL));

    printf("========================================================\n");
    printf(" 启动固件级 LDPC 引擎能力探测仿真\n");
    printf(" 配置: Data=%d Bytes, Parity=%d Bytes (Rate: %.2f)\n", DATA_BYTES,
           PARITY_BYTES, (float)DATA_BYTES / (DATA_BYTES + PARITY_BYTES));
    printf("========================================================\n");

    /* 初始化 */
    #if 0
    build_mock_rom_table();
    if (ldpc_engine_init() < 0)
    {
        printf("[Error] LDPC 引擎初始化失败！\n");
        return -1;
    }
    #elif 0
    /* 召唤寻优器！让它穷举 50 种矩阵结构，选出最好的一套 */
    find_optimal_matrix(50000);
    #else
    build_mock_rom_table_qc();
    if (ldpc_engine_init() < 0)
    {
        printf("[Error] LDPC 引擎初始化失败！\n");
        return -1;
    }
    #endif

    /* 随机生成宿主机数据 */
    srand((unsigned)time(NULL));
    for (i = 0; i < DATA_BITS; i++)
    {
        raw_data[i] = rand() % 2;
    }

    /* 编码落盘 */
    printf("\n[FW Process] 正在执行系统码编码...\n");
    ldpc_encode(raw_data, tx_codeword);

    /* 运行压力测试 */
    test_correction_limit(tx_codeword);

    /* --- 新增：运行宏观系统容灾演习 --- */
    firmware_raid_recovery_demo();

    /* 清理退场 */
    printf("\n[FW Process] 仿真完毕，释放资源。\n");
    free(raw_data);
    free(tx_codeword);
    ldpc_engine_cleanup();

    return 0;
}