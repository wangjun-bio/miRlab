import argparse
import gzip
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

# --------------------------
# 函数：每块处理
# --------------------------
def process_chunk(chunk):
    result = []
    for record in chunk:
        result.append(f"{record.id}\t{len(record.seq)}")
    return result

# --------------------------
# 主程序入口
# --------------------------
def main():
    parser = argparse.ArgumentParser(description="并行提取 FASTQ 序列长度")
    parser.add_argument("-i", "--input", required=True, help="输入 FASTQ(.gz) 文件")
    parser.add_argument("-o", "--output", required=True, help="输出 TSV 文件")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count(), help="使用的线程数（默认最大）")
    parser.add_argument("-b", "--batch", type=int, default=10000, help="每个进程处理的记录数（默认10000）")
    args = parser.parse_args()

    # 自动选择打开方式
    open_func = gzip.open if args.input.endswith(".gz") else open

    # 分块读取
    with open_func(args.input, "rt") as handle, Pool(args.threads) as pool, open(args.output, "w") as out:
        batch = []
        for record in SeqIO.parse(handle, "fastq"):
            batch.append(record)
            if len(batch) >= args.batch:
                for line in pool.map(process_chunk, [batch]):
                    out.write("\n".join(line) + "\n")
                batch = []

        # 处理最后不足一块的数据
        if batch:
            for line in pool.map(process_chunk, [batch]):
                out.write("\n".join(line) + "\n")

if __name__ == "__main__":
    main()
