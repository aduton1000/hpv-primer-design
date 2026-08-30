import shutil
import os
#!/usr/bin/env python
"""
v4 TASK 6 - re-confirm the human (GRCh38) screen is unchanged: 0 predicted human
amplicons. Only HPV58_3 primers changed vs v3; all others are identical. Same logic
as v3 (both primers 3'-anchored, opposite strands, convergent, <=5 kb on the same
human chromosome -> linear distance, human genome is not circular).
"""
import subprocess
from pathlib import Path
from collections import defaultdict

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
BL = os.environ.get("HPV_BLAST_BIN") or (os.path.dirname(shutil.which("blastn")) if shutil.which("blastn") else "")
GRCH38 = PROJ / "v3" / "01_refs" / "grch38"      # full GRCh38 toplevel (706 seqs, 3.29 Gb)
DES = PROJ / "v5" / "04_design"
SPEC = PROJ / "v5" / "05_specificity"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
MAX_AMP = 5000
log = []
def L(m): print(m); log.append(str(m))


def load_primers():
    prim = {}; amps = defaultdict(lambda: defaultdict(dict))
    for t in TYPES:
        for ln in (DES / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t")
            prim[c[3]] = dict(type=f"HPV{t}", amp=c[1], dir=c[4], seq=c[5], length=len(c[5]))
            amps[f"HPV{t}"][c[1]][c[4]] = c[3]
    return prim, amps


def main():
    prim, amps = load_primers()
    plen = {n: p["length"] for n, p in prim.items()}
    pf = SPEC / "all_primers.fasta"
    pf.write_text("\n".join(f">{n}\n{p['seq']}" for n, p in prim.items()) + "\n")
    L(f"=== T6 human screen (GRCh38, 706 seqs / 3.29 Gb), {len(prim)} primers ===")
    fmt = ("6 qseqid sseqid pident length qlen qstart qend sstart send sstrand "
           "mismatch gapopen qseq sseq")
    r = subprocess.run([f"{BL}/blastn", "-query", str(pf), "-db", str(GRCH38),
                        "-task", "blastn-short", "-evalue", "100", "-word_size", "7",
                        "-outfmt", fmt, "-max_target_seqs", "50000", "-num_threads", "6"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        L(f"blast rc={r.returncode}: {r.stderr[:300]}"); return None
    hits = []
    for ln in r.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 14: continue
        q = f[0]; pl = plen.get(q, int(f[4]))
        length = int(f[3]); qend = int(f[6]); mism = int(f[10]) + int(f[11])
        qseq, sseq = f[12], f[13]
        tp = (qend == pl and len(qseq) >= 5 and qseq[-5:] == sseq[-5:])
        binds = (length >= pl - 3 and mism <= 3 and tp)
        hits.append(dict(primer=q, subj=f[1], strand=f[9], sstart=int(f[7]), send=int(f[8]),
                         binds=binds, pident=float(f[2]), mism=mism, qend=qend, qlen=pl))
    binders = sorted({h["primer"] for h in hits if h["binds"]})
    with open(SPEC / "human_hits.tsv", "w") as fh:
        fh.write("primer\tchrom\tpident\tmismatch\tqend\tqlen\tstrand\tsstart\tsend\t3prime_anchored_bind\n")
        for h in hits:
            if h["binds"]:
                fh.write(f"{h['primer']}\t{h['subj']}\t{h['pident']}\t{h['mism']}\t{h['qend']}\t"
                         f"{h['qlen']}\t{h['strand']}\t{h['sstart']}\t{h['send']}\tYES\n")
    # amplicon caller (both primers, opposite strand, convergent, <=5kb on same chrom)
    by_subj = defaultdict(list)
    for h in hits:
        if h["binds"]:
            by_subj[h["subj"]].append(h)
    amp_rows = []
    for typ, ad in amps.items():
        for amp, fr in ad.items():
            fn, rn = fr.get("F"), fr.get("R")
            for subj, sh in by_subj.items():
                fs = [h for h in sh if h["primer"] == fn]
                rs = [h for h in sh if h["primer"] == rn]
                for a in fs:
                    for b in rs:
                        if a["strand"] == b["strand"]: continue
                        plus = a if a["strand"] == "plus" else b
                        minus = b if a["strand"] == "plus" else a
                        gap = minus["send"] - plus["send"]
                        if 50 <= gap <= MAX_AMP:
                            amp_rows.append((typ, amp, subj, gap))
    with open(SPEC / "human_amplicons.tsv", "w") as fh:
        fh.write("type\tamplicon\tchrom\tpredicted_product_bp\n")
        for r2 in amp_rows:
            fh.write("\t".join(str(x) for x in r2) + "\n")
    L(f"  human 3'-anchored single-primer binders: {len(binders)}")
    L(f"  human predicted AMPLICONS (both primers, convergent, <=5kb): {len(amp_rows)}")
    # did HPV58_3 (the only changed primers) introduce anything?
    h58 = [b for b in binders if b.startswith("HPV58_3_")]
    L(f"  HPV58_3 (redesigned) human binders: {h58 if h58 else 'none'}")
    return len(binders), len(amp_rows)


if __name__ == "__main__":
    nb, na = main()
    (PROJ / "v5" / "00_logs" / "t5_human.log").write_text("\n".join(log) + "\n")
    print(f"Log -> {PROJ/'v4/00_logs/t6_human.log'}")
