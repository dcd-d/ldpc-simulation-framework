import gdb

class PlotLDPCSquare(gdb.Command):
    """打印单个 QC-LDPC 宏块展开后的比特矩阵"""
    def __init__(self):
        super(PlotLDPCSquare, self).__init__("plot_block", gdb.COMMAND_USER)

    def invoke(self, arg, from_tty):
        args = gdb.string_to_argv(arg)
        shift = int(gdb.parse_and_eval(args[0]))
        z = int(gdb.parse_and_eval(args[1]))

        print(f"\n[Matrix Block] Shift={shift}, Z={z}")
        header = "   " + " ".join([f"{j%10}" for j in range(z)])
        print(header)
        print("  " + "-" * (len(header)-2))

        for i in range(z):
            target_j = (i + shift) % z
            row = f"{i%10}| " + " ".join(["1" if j == target_j else "." for j in range(z)])
            print(row)

class PlotCodeword(gdb.Command):
    """打印整个 LDPC 码字（数据位 + 校验位）"""
    def __init__(self):
        super(PlotCodeword, self).__init__("plot_codeword", gdb.COMMAND_USER)

    def invoke(self, arg, from_tty):
        args = gdb.string_to_argv(arg)
        ptr = gdb.parse_and_eval(args[0])
        total_bits = int(gdb.parse_and_eval(args[1]))
        data_bits = 64 # 根据你的迷你配置

        print(f"\n[Codeword Visualization] Total: {total_bits} bits")
        print("Data Area " + "-" * 30 + " | Parity Area")
        
        output = ""
        for i in range(total_bits):
            if i == data_bits: output += " | " # 分隔符
            bit = int(ptr[i])
            output += "1" if bit else "0"
            if (i+1) % 8 == 0 and i+1 != data_bits and i+1 != total_bits:
                output += " " # 字节分隔
        
        print(output + "\n")

PlotLDPCSquare()
PlotCodeword()