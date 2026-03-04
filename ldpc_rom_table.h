#ifndef LDPC_ROM_TABLE_H
#define LDPC_ROM_TABLE_H

/* 物理参数定义 */
#define DATA_BYTES              4096
#define PARITY_BYTES            256
#define DATA_BITS               (DATA_BYTES * 8)
#define PARITY_BITS             (PARITY_BYTES * 8)
#define CODEWORD_BITS           (DATA_BITS + PARITY_BITS)

#define LIFTING_FACTOR          64
#define NUM_DATA_BLOCKS         (DATA_BITS / LIFTING_FACTOR)
#define NUM_CHK_BLOCKS          (PARITY_BITS / LIFTING_FACTOR)

#define MAX_ITERATIONS          20
#define MAX_HARD_ITERS          100

/* * 稀疏矩阵坐标条目结构
 * 在实际芯片中，这通常对应于内部连线路由表的配置寄存器
 */
struct ldpc_rom_entry {
        unsigned short row;
        unsigned short col;
        unsigned short shift;
};

/* * 模拟算法团队交付的最优 QC-LDPC 基矩阵参数表
 * 注意：为了演示，这里仅放置了少量示例数据。
 * 在真实工程中，这个数组可能有数千行，由算法团队生成的脚本自动导出。
 */
static const struct ldpc_rom_entry g_rom_base_matrix[] = {
        {0, 0, 15}, {0, 15, 32}, {0, 256, 7},  {0, 500, 61}, 
        {1, 1, 44}, {1, 16, 12}, {1, 257, 55}, {1, 501, 3},  
        /* ... 算法团队提供的成百上千个最优坐标点 ... */
        {31, 14, 8}, {31, 255, 19}, {31, 499, 42}, {31, 511, 0}
};

/* 计算表项总数 */
#define ROM_MATRIX_ENTRIES (sizeof(g_rom_base_matrix) / sizeof(g_rom_base_matrix[0]))

#endif /* LDPC_ROM_TABLE_H */