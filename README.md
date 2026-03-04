# Enterprise SSD LDPC & RAIN Simulation Framework
**企业级固件 LDPC 纠错与跨 Die 容灾仿真引擎 (C 语言实现)**



## 📖 项目简介 (Overview)
本项目是一个基于 C 语言构建的底层固件级闪存纠错与数据保护仿真框架。它不仅实现了标准的 LDPC (低密度奇偶校验) 编解码算法，还深度还原了真实 3D NAND 闪存的物理特性与企业级 SSD 主控芯片的异常处理流程。

从物理层的 **3-bit 软读取量化**，到算法层的 **无 4 环 QC-LDPC 矩阵生成**，再到系统架构层的 **跨 Die 异或容灾 (RAIN)**，全链路展示了现代固件保护用户数据的极致防御体系。

## ✨ 核心特性 (Key Features)

### 1. 算法级图模型优化：QC-LDPC 无短环生成器
抛弃传统的随机连线，采用代数同余方程构造准循环 (Quasi-Cyclic) 矩阵。
通过严格的数学约束彻底消灭图模型中的短环 (Cycle-4)：
$S_{11} + S_{22} \equiv S_{12} + S_{21} \pmod Z$
* **效果**：极大地拓宽了硬判决 (Hard Decode) 的纠错边界，消除了早期 Error Floor。

### 2. 真实物理层建模：3-bit NAND Soft Read 查表法 (LUT)
抛弃了理想化的双精度浮点 LLR 计算，引入了符合真实硬件设计的 3-bit 阈值电压切片模型。
* **效果**：通过位运算拼凑离散的置信度索引，精准还原了硬件量化带来的精度截断误差 (Quantization Loss)，在仿真中复现了真实的单 Page 纠错物理极限。

### 3. 防死锁硬判决：概率比特翻转 (PBF / GDBF)
在贪心比特翻转 (Max-BF) 的基础上，引入了宏观伴随式权重监控 (Syndrome Weight Tracking) 与死锁状态机。
* **效果**：当检测到解码器陷入陷阱集 (Trapping Sets) 时，引入伪随机噪声进行概率翻转，利用局部混沌打破对称性死锁。

### 4. 工业级软判决：归一化最小和算法 (NMSA)
实现了企业级主控标准的 Normalized Min-Sum 算法，包含 0.75 衰减系数与提前退出机制 (Early Termination)，在纠错能力与执行功耗之间取得完美平衡。

### 5. 终极数据防线：跨 Die 异或容灾重构 (RAIN)
模拟了真实 SSD 固件的极限求生状态机：
当某一个物理 Die 遭遇致死级漏电 (如 SNR 跌至 2.5dB 触发 LDPC UE) 时，固件自动挂起 Host I/O，并发读取冗余 Die (RAID 5 架构)，通过底层 XOR 重构完美恢复用户数据，避免蓝屏 (BSOD)。

## 🚀 快速开始 (Getting Started)

本项目无任何第三方依赖，纯标准 C 语言编写，支持在 Linux (GCC/Clang) 或 Windows (MinGW) 下一键编译。

### 编译与运行
```bash
# 1. 克隆代码库
git clone [https://github.com/](https://github.com/)<Your-Username>/ldpc-simulation-framework.git
cd ldpc-simulation-framework

# 2. 编译并开启 GDB 调试支持
gcc -o ldpc_simu ldpc_decoder.c -lm -O0 -ggdb

# 3. 启动仿真引擎
./ldpc_simu