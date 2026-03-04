// #define DATA_BYTES 8
#define DATA_BYTES 4096

#if DATA_BYTES == 4096
#define DATA_BYTES 4096
#define PARITY_BYTES 256
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
#define DATA_BYTES              8
#define PARITY_BYTES            2

#define DATA_BITS               (DATA_BYTES * 8)             // 64 bits
#define PARITY_BITS             (PARITY_BYTES * 8)           // 16 bits
#define CODEWORD_BITS           (DATA_BITS + PARITY_BITS)    // 80 bits

#define LIFTING_FACTOR          4                            // 矩阵由 4x4 的小块组成
#define NUM_DATA_BLOCKS         (DATA_BITS / LIFTING_FACTOR) // 16 列
#define NUM_CHK_BLOCKS          (PARITY_BITS / LIFTING_FACTOR)// 4 行
#define COLUMN_WEIGHT           3                            // 每列 3 个非零块

#define MAX_ITERATIONS          100
#define MAX_HARD_ITERS          5000
#endif