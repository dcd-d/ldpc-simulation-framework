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

    srand(12345); // 固定随机种子，保证每次调试矩阵结构一致
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
    printf("[Testbed] 成功生成模拟 ROM 表，包含 %d 个有效宏块。\n",
           g_rom_matrix_entries);
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
    printf("[LDPC Engine] 硬件图模型初始化完毕。\n");
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
/* ==========================================
 * 模块 3: 核心解码算法
 * ========================================== */
int ldpc_decode_hard(unsigned char *rx_codeword)
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
 * 主函数入口
 * ========================================== */
int main(void)
{
    unsigned char *raw_data = (unsigned char *)malloc(DATA_BITS);
    unsigned char *tx_codeword = (unsigned char *)malloc(CODEWORD_BITS);
    int i;

    printf("========================================================\n");
    printf(" 启动固件级 LDPC 引擎能力探测仿真\n");
    printf(" 配置: Data=%d Bytes, Parity=%d Bytes (Rate: %.2f)\n", DATA_BYTES,
           PARITY_BYTES, (float)DATA_BYTES / (DATA_BYTES + PARITY_BYTES));
    printf("========================================================\n");

    /* 初始化 */
    build_mock_rom_table();
    if (ldpc_engine_init() < 0)
    {
        printf("[Error] LDPC 引擎初始化失败！\n");
        return -1;
    }

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

    /* 清理退场 */
    printf("\n[FW Process] 仿真完毕，释放资源。\n");
    free(raw_data);
    free(tx_codeword);
    ldpc_engine_cleanup();

    return 0;
}