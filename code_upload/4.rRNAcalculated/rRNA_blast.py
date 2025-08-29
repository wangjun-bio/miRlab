#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path
from subprocess import Popen, PIPE, CalledProcessError
import glob
import shutil
import tempfile
from typing import List, Tuple
from tqdm import tqdm


# ========== 工具函数 ==========

def check_blast_db(prefix: str) -> bool:
    """判断 BLAST 数据库前缀是否存在（常见索引文件之一存在即可）"""
    patterns = [prefix + ext for ext in [".nhr", ".nin", ".nsq", ".nal", ".phr", ".pin", ".psq"]]
    return any(glob.glob(p) for p in patterns)


def run_and_check(cmd, **popen_kwargs):
    """运行子进程并校验返回码，出错抛出 CalledProcessError"""
    proc = Popen(cmd, **popen_kwargs)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise CalledProcessError(proc.returncode, cmd, output=stdout, stderr=stderr)
    return stdout, stderr


def sample_one_bam_to_fasta(
    bam: Path,
    out_fa: Path,
    samtools_threads: int,
    per_bam: int,
    seed: int
):
    """
    从单个 BAM 抽样未比对 reads 的 N 条，转为 FASTA 写入 out_fa
    命令链：samtools fastq -f 4 | seqtk sample -s SEED N | seqtk seq -A -
    """
    # 1) samtools 未比对 reads → FASTQ
    p1 = Popen(
        ["samtools", "fastq", "-@", str(samtools_threads), "-f", "4", str(bam)],
        stdout=PIPE, stderr=PIPE, text=False
    )
    # 2) 随机抽样 N 条
    p2 = Popen(
        ["seqtk", "sample", "-s", str(seed), "-", str(per_bam)],
        stdin=p1.stdout, stdout=PIPE, stderr=PIPE, text=False
    )
    if p1.stdout: p1.stdout.close()
    # 3) 转为 FASTA
    p3 = Popen(
        ["seqtk", "seq", "-A", "-"],
        stdin=p2.stdout, stdout=PIPE, stderr=PIPE, text=False
    )
    if p2.stdout: p2.stdout.close()

    with open(out_fa, "wb") as fout:
        data, _ = p3.communicate()
        if data:
            fout.write(data)

    # 等待并检查前级
    _, err3 = p3.communicate() if p3.poll() is None else (b"", b"")
    rc3 = p3.returncode
    _, err2 = p2.communicate()
    rc2 = p2.returncode
    _, err1 = p1.communicate()
    rc1 = p1.returncode

    if rc1 != 0:
        raise CalledProcessError(rc1, "samtools fastq", stderr=err1)
    if rc2 != 0:
        raise CalledProcessError(rc2, "seqtk sample", stderr=err2)
    if rc3 != 0:
        raise CalledProcessError(rc3, "seqtk seq -A", stderr=err3)


def count_fasta_seqs(fa_path: Path) -> int:
    """统计 FASTA 文件中的序列条数（以 '>' 开头的行数）"""
    n = 0
    with open(fa_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def split_fasta_to_chunks(
    fa_path: Path,
    out_dir: Path,
    chunk_size: int
) -> List[Tuple[Path, int]]:
    """
    按序列条数把 FASTA 切分成多个小文件（每块最多 chunk_size 条）
    返回 [(chunk_path, count), ...]
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks: List[Tuple[Path, int]] = []
    idx = 0
    cur_count = 0
    cur_path = out_dir / f"chunk_{idx:05d}.fa"
    cur_out = open(cur_path, "wt")
    chunks.append((cur_path, 0))  # 先占位，稍后更新 count

    def _new_chunk():
        nonlocal idx, cur_count, cur_path, cur_out
        cur_out.close()
        # 更新上一块的 count
        chunks[-1] = (cur_path, cur_count)
        # 开启新块
        idx += 1
        cur_count = 0
        cur_path = out_dir / f"chunk_{idx:05d}.fa"
        cur_out = open(cur_path, "wt")
        chunks.append((cur_path, 0))

    with open(fa_path, "rt") as f:
        writing = False
        for line in f:
            if line.startswith(">"):
                # 新序列开始：如当前块已满，换新块
                if cur_count >= chunk_size:
                    _new_chunk()
                cur_count += 1
                writing = True
            if writing:
                cur_out.write(line)

    # 关闭最后一块并更新计数
    cur_out.close()
    chunks[-1] = (cur_path, cur_count)

    # 去掉可能的空块
    chunks = [(p, c) for (p, c) in chunks if c > 0]
    return chunks


def blast_chunk_outfmt6(
    query_fa: Path,
    db_prefix: str,
    out_handle,
    task: str,
    evalue: str,
    perc_identity: str,
    qcov_hsp_perc: str,
    max_hsps: str,
    threads: int,
):
    """
    对单个 FASTA 块运行 blastn（outfmt 6），结果写入已经打开的 out_handle（以追加方式）
    """
    outfmt = "6 qseqid sacc evalue bitscore pident stitle"
    cmd = [
        "blastn",
        "-task", task,
        "-query", str(query_fa),
        "-db", db_prefix,                 # 前缀作为字符串直接传入
        "-evalue", str(evalue),
        "-perc_identity", str(perc_identity),
        "-qcov_hsp_perc", str(qcov_hsp_perc),
        "-max_hsps", str(max_hsps),
        "-num_threads", str(threads),
        "-outfmt", outfmt,
    ]
    # 直接把 stdout 追加到 out_handle
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, text=True, bufsize=1)
    for line in proc.stdout:
        out_handle.write(line)
    stderr = proc.stderr.read()
    ret = proc.wait()
    if ret != 0:
        raise CalledProcessError(ret, cmd, stderr=stderr)


# ========== 主流程 ==========

def main():
    ap = argparse.ArgumentParser(
        description="从每个 BAM 抽样 N 条未比对 reads，合并后按块运行 BLASTN（outfmt 6）并显示准确进度"
    )
    ap.add_argument("--input-dir", default="./SILVA_NR99_STAR_LSU", help="BAM 所在目录")
    ap.add_argument("--glob", default="*_Aligned.sortedByCoord.out.bam", help="BAM 匹配模式")
    ap.add_argument("--per-bam", type=int, default=100, help="每个 BAM 抽取的 reads 数 (默认 100)")
    ap.add_argument("--seed", type=int, default=42, help="随机种子（用于 seqtk sample）")
    ap.add_argument("--samtools-threads", type=int, default=4, help="samtools 线程数")
    ap.add_argument("--db", required=True, help="BLAST 数据库前缀（字符串原样传入，不当作目录）")
    ap.add_argument("--outdir", default="blastn_out", help="输出目录")
    ap.add_argument("--outfile", default="combined_sampled.out", help="最终 BLAST 输出文件名（outfmt 6）")
    ap.add_argument("--task", default="megablast", choices=["megablast","blastn","dc-megablast"], help="blastn -task")
    ap.add_argument("--evalue", default="1e-5", help="E-value 阈值")
    ap.add_argument("--perc-identity", default="0", help="最小百分比相似度")
    ap.add_argument("--qcov-hsp-perc", default="0", help="query 覆盖率阈值（HSP 级）")
    ap.add_argument("--max-hsps", default="1", help="每个 subject 的 HSP 数")
    ap.add_argument("--blastn-threads", type=int, default=8, help="blastn 线程数")
    ap.add_argument("--chunk-size", type=int, default=1000, help="合并 FASTA 后每块包含的序列条数（用于进度条）")
    ap.add_argument("--keep-temp", action="store_true", help="保留中间临时文件（抽样块、合并 FA、切块 FA）")
    ap.add_argument("--tmpdir", default=None, help="指定临时文件目录（默认使用系统临时目录）")
    args = ap.parse_args()

    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) 收集 BAM
    bam_files = sorted(input_dir.glob(args.glob))
    if not bam_files:
        print(f"[ERR] 未找到 BAM：{args.glob} in {input_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] 发现 BAM 数量：{len(bam_files)}")

    # 2) 检查 BLAST DB（保持前缀字符串，不当作目录）
    db_prefix = args.db
    if not check_blast_db(db_prefix):
        print(f"[ERR] 未找到 BLAST DB 前缀对应文件：{db_prefix}.(nhr/nin/nsq/…)", file=sys.stderr)
        sys.exit(2)

    # 3) 抽样每个 BAM → 临时抽样 FASTA
    if args.tmpdir:
        tmp_root_obj = tempfile.TemporaryDirectory(dir=args.tmpdir)
    else:
        tmp_root_obj = tempfile.TemporaryDirectory()
    tmp_root = Path(tmp_root_obj.name)

    tmp_chunks_dir = tmp_root / "chunks"
    tmp_samples_dir = tmp_root / "samples"
    tmp_samples_dir.mkdir(parents=True, exist_ok=True)

    sample_paths = []
    for i, bam in enumerate(bam_files, 1):
        out_fa = tmp_samples_dir / f"{bam.stem}.sample{args.per_bam}.fa"
        print(f"[INFO] 抽样 ({i}/{len(bam_files)}): {bam.name} -> {out_fa.name}")
        try:
            sample_one_bam_to_fasta(
                bam=bam,
                out_fa=out_fa,
                samtools_threads=args.samtools_threads,
                per_bam=args.per_bam,
                seed=args.seed
            )
            sample_paths.append(out_fa)
        except CalledProcessError as e:
            msg = (e.stderr.decode("utf-8", errors="ignore")
                   if isinstance(e.stderr, (bytes, bytearray)) else str(e))
            print(f"[WARN] 抽样失败，跳过 {bam.name}: {msg}", file=sys.stderr)

    if not sample_paths:
        print("[ERR] 没有可用的抽样序列，终止。", file=sys.stderr)
        sys.exit(3)

    # 4) 合并为一个大 FASTA
    combined_fa = outdir / f"combined_{args.per_bam}perbam.fa"
    print(f"[INFO] 合并 {len(sample_paths)} 个抽样 FASTA -> {combined_fa}")
    with open(combined_fa, "wb") as fout:
        for p in sample_paths:
            with open(p, "rb") as fin:
                shutil.copyfileobj(fin, fout)

    total_queries = count_fasta_seqs(combined_fa)
    if total_queries == 0:
        print("[ERR] 合并后 FASTA 序列数为 0。", file=sys.stderr)
        sys.exit(4)
    print(f"[INFO] 合并后 query 总数：{total_queries}")

    # 5) 切块（每块 chunk_size 条）
    chunks = split_fasta_to_chunks(combined_fa, tmp_chunks_dir, args.chunk_size)
    print(f"[INFO] 生成 {len(chunks)} 个 chunk（每块最多 {args.chunk_size} 条）")

    # 6) 逐块运行 blastn（outfmt 6），并用进度条累计精确的已处理 query 数
    out_path = outdir / args.outfile
    # 先清空/创建输出文件
    open(out_path, "wt").close()

    processed = 0
    with open(out_path, "a") as fout, tqdm(total=total_queries, desc="blastn", unit="query") as pbar:
        for chunk_path, n_in_chunk in chunks:
            try:
                blast_chunk_outfmt6(
                    query_fa=chunk_path,
                    db_prefix=db_prefix,
                    out_handle=fout,
                    task=args.task,
                    evalue=args.evalue,
                    perc_identity=args.perc_identity,
                    qcov_hsp_perc=args.qcov_hsp_perc,
                    max_hsps=args.max_hsps,
                    threads=args.blastn_threads
                )
                processed += n_in_chunk
                pbar.update(n_in_chunk)
            except CalledProcessError as e:
                msg = e.stderr if isinstance(e.stderr, str) else str(e)
                print(f"[ERR] BLAST 失败（{chunk_path.name}）: {msg}", file=sys.stderr)
                sys.exit(5)

    print(f"[DONE] 完成。结果：{out_path} ；总计处理 {processed}/{total_queries} 条 query")

    # 7) 是否保留临时文件
    if args.keep_temp:
        print(f"[KEEP] 临时文件保留于：{tmp_root}")
    else:
        tmp_root_obj.cleanup()


if __name__ == "__main__":
    main()
