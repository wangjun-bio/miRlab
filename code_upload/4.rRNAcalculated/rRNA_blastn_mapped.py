import os, re, gzip
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ----- 关键词定义（精确到词边界）；顺序：先具体，最后泛化的 'rRNA' -----
TERMS = [
    "5S", "16S", "18S", "23S", "28S",
    "SSU rRNA", "LSU rRNA", "ribosomal RNA", "rRNA" ,'5.8S'
]
# 为每个关键词编译独立正则（大小写不敏感）
# 使用 \b 词边界（等价于 (?<!\w) 与 (?!\w) 的简写），避免部分词误匹配
TERM_PATS = [(t, re.compile(rf"\b{re.escape(t)}\b", re.I)) for t in TERMS]

def iter_lines(path: Path):
    if str(path).endswith(".gz"):
        f = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        f = open(path, "rt", encoding="utf-8", newline="", buffering=1024 * 1024)
    with f:
        for line in f:
            yield line.rstrip("\n")

def match_terms(text: str) -> set:
    """
    返回该行命中的关键词集合；若命中更具体的 rRNA 词条，则去掉泛化的 'rRNA'。
    """
    hits = {t for t, pat in TERM_PATS if pat.search(text)}
    # 若已命中更具体的 rRNA 描述，则去掉泛化的 'rRNA'，减少冗余
    specific = {"SSU rRNA", "LSU rRNA", "ribosomal RNA"}
    if hits & specific and "rRNA" in hits:
        hits.remove("rRNA")
    return hits

def summarize_file(path: Path):
    """
    逐文件处理：
      - 建立 qseqid -> set(命中关键词) 的映射
      - 返回 (文件名，命中数，未命中数，细节列表 [ (qseqid, 'term1;term2') ... ])
    """
    q2terms = {}  # qseqid -> set of terms
    for line in iter_lines(path):
        if not line:
            continue
        parts = line.split("\t")
        if not parts:
            continue
        q = parts[0]
        # 初始化集合
        s = q2terms.get(q)
        if s is None:
            s = set()
            q2terms[q] = s
        # 抽取该行命中的关键词，合并到该 qseqid 的集合
        s |= match_terms(line)

    # 统计匹配/未匹配
    matched = sum(1 for terms in q2terms.values() if terms)
    unmatched = sum(1 for terms in q2terms.values() if not terms)

    # 细节行：仅输出命中的 qseqid；如需包含未命中，可在此分支中加入 NA
    details = []
    for q, terms in q2terms.items():
        if terms:
            details.append((q, ";".join(sorted(terms))))
    return path.name, matched, unmatched, details

def main(src, pattern="*.out", summary="summary.tsv", details="details.tsv", workers=8):
    files = [p for p in Path(src).rglob(pattern) if p.is_file() and p.stat().st_size > 0]
    with ProcessPoolExecutor(max_workers=workers) as ex, \
         open(summary, "w") as fsum, \
         open(details, "w") as fdet:
        # 写表头
        fsum.write("file\tmatched_rrna_qseqids\tunmatched_qseqids\n")
        fdet.write("file\tqseqid\tmatched_terms\n")

        futures = {ex.submit(summarize_file, p): p for p in files}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
            fn, matched, unmatched, det_rows = fut.result()
            # 写 summary
            fsum.write(f"{fn}\t{matched}\t{unmatched}\n")
            # 写 details（逐 qseqid 列出命中的关键词）
            if det_rows:
                for q, terms in det_rows:
                    fdet.write(f"{fn}\t{q}\t{terms}\n")

if __name__ == "__main__":
    main(
        src="/mnt/4T_SSD/yiyonghao/blastn_out",
        pattern="*.out",
        summary="summary.tsv",
        details="details.tsv",
        workers=30
    )
