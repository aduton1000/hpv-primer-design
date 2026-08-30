import shutil
import os
#!/usr/bin/env python
"""
v4 TASK 3 - Off-target (cross-organism) amplicon-level in-silico PCR.

Panel = v3's extra panel (11 non-panel HPV + 12 cervicovaginal/other organisms)
PLUS Trichomonas vaginalis G3 (GCF_000002825.2 / GCA_000002825.1, ASM282v1,
assembly level=Scaffold, WGS AAHC01, 64,769 scaffolds, 176.4 Mb) -- the organism
that was UNAVAILABLE as a complete RefSeq in v3 and is now included as a draft.

Same amplicon logic as the human screen: a hit counts only if BOTH primers of an
owner amplicon 3'-anchor (3' terminal 5 nt EXACT, <=3 total mismatch, full-length)
on the SAME subject sequence, on OPPOSITE strands, in convergent orientation,
within a plausible product distance (<=5 kb). Off-target organism genomes are linear
(or independent scaffolds), so distance is linear per-subject.

Outputs (v4/05_specificity/):
  offtarget_hits.tsv      primer, owner_type, hit_subject, organism, sstart, send,
                          strand, pident, mismatch, 3p_bind
  offtarget_amplicons.tsv owner_amplicon, organism, hit_subject, predicted_product_bp, verdict
"""
import subprocess, sys
from pathlib import Path
from collections import defaultdict

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
BL = os.environ.get("HPV_BLAST_BIN") or (os.path.dirname(shutil.which("blastn")) if shutil.which("blastn") else "")
DB = PROJ / "v4" / "01_refs" / "extra"
SPEC = PROJ / "v5" / "05_specificity"
DES = PROJ / "v5" / "04_design"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
MAX_AMP = 5000
log = []
def L(m): print(m); log.append(str(m))


def load_primers(pdir):
    prim = {}; amps = defaultdict(lambda: defaultdict(dict))
    for t in TYPES:
        for ln in (Path(pdir) / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t")
            prim[c[3]] = dict(type=f"HPV{t}", amp=c[1], dir=c[4], seq=c[5], length=len(c[5]))
            amps[f"HPV{t}"][c[1]][c[4]] = c[3]
    return prim, amps


def blast(query_fa, db, plen):
    fmt = ("6 qseqid sseqid pident length qlen qstart qend sstart send sstrand "
           "mismatch gapopen qseq sseq")
    r = subprocess.run([f"{BL}/blastn", "-query", str(query_fa), "-db", str(db),
                        "-task", "blastn-short", "-evalue", "100", "-word_size", "7",
                        "-outfmt", fmt, "-max_target_seqs", "100000", "-num_threads", "6"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        L(f"  blast rc={r.returncode}: {r.stderr[:300]}"); return []
    hits = []
    for ln in r.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 14: continue
        q = f[0]; pl = plen.get(q, int(f[4]))
        length = int(f[3]); qend = int(f[6]); mism = int(f[10]) + int(f[11])
        qseq, sseq = f[12], f[13]
        three_prime_exact = (qend == pl and len(qseq) >= 5 and qseq[-5:] == sseq[-5:])
        binds = (length >= pl - 3 and mism <= 3 and three_prime_exact)
        hits.append(dict(primer=q, subj=f[1], strand=f[9], sstart=int(f[7]), send=int(f[8]),
                         binds=binds, pident=float(f[2]), length=length, mism=mism,
                         qend=qend, qlen=pl))
    return hits


def main():
    pdir = sys.argv[1] if len(sys.argv) > 1 else str(DES)
    prim, amps = load_primers(pdir)
    plen = {n: p["length"] for n, p in prim.items()}
    pf = SPEC / "all_primers.fasta"
    pf.write_text("\n".join(f">{n}\n{p['seq']}" for n, p in prim.items()) + "\n")
    L(f"=== T3 off-target screen ({len(prim)} primers) vs combined panel ===")
    L(f"  db = {DB} (35 v3 panel seqs + 64,769 T. vaginalis scaffolds = 64,804)")

    hits = blast(pf, DB, plen)
    L(f"  total raw BLAST hits: {len(hits)}")

    # hits table (binds + near-misses for audit)
    nbind = 0
    org_binds = defaultdict(int)
    with open(SPEC / "offtarget_hits.tsv", "w") as fh:
        fh.write("primer\towner_type\thit_subject\torganism\tsstart\tsend\tstrand\tpident\t"
                 "mismatch\t3p_bind\n")
        for h in sorted(hits, key=lambda x: (x["primer"], x["subj"])):
            org = h["subj"].split("__")[0]
            if h["binds"] or (h["length"] >= h["qlen"] - 3 and h["mism"] <= 4):
                fh.write(f"{h['primer']}\t{prim[h['primer']]['type']}\t{h['subj']}\t{org}\t"
                         f"{h['sstart']}\t{h['send']}\t{h['strand']}\t{h['pident']}\t{h['mism']}\t"
                         f"{'YES' if h['binds'] else 'no'}\n")
                if h["binds"]:
                    nbind += 1; org_binds[org] += 1
    L(f"  3'-anchored single-primer binds: {nbind}")
    if org_binds:
        L(f"    by organism: {dict(sorted(org_binds.items()))}")

    # amplicon caller (both primers, same subject, opposite strand, convergent, <=5kb)
    by_subj = defaultdict(list)
    for h in hits:
        if h["binds"]:
            by_subj[h["subj"]].append(h)
    amp_rows = []
    for typ, ad in amps.items():
        for amp, fr in ad.items():
            fn, rn = fr.get("F"), fr.get("R")
            if not fn or not rn: continue
            for subj, sh in by_subj.items():
                fs = [h for h in sh if h["primer"] == fn]
                rs = [h for h in sh if h["primer"] == rn]
                for a in fs:
                    for b in rs:
                        if a["strand"] == b["strand"]:
                            continue
                        plus = a if a["strand"] == "plus" else b
                        minus = b if a["strand"] == "plus" else a
                        gap = minus["send"] - plus["send"]   # linear (off-target genomes/scaffolds)
                        if 50 <= gap <= MAX_AMP:
                            amp_rows.append((f"{typ}_{amp}", subj.split("__")[0], subj, gap))

    with open(SPEC / "offtarget_amplicons.tsv", "w") as fh:
        fh.write("owner_amplicon\torganism\thit_subject\tpredicted_product_bp\tverdict\n")
        for r in amp_rows:
            fh.write("\t".join(str(x) for x in r) + "\tFAIL\n")
    L(f"  OFF-TARGET AMPLICONS (both-primer, convergent, <=5kb): {len(amp_rows)}")
    for r in amp_rows:
        L(f"    FAIL: {r[0]} -> {r[1]} ({r[2]}) {r[3]} bp")

    # explicit T. vaginalis verdict
    tvag_binds = org_binds.get("Trichomonas_vaginalis", 0)
    tvag_amps = sum(1 for r in amp_rows if r[1] == "Trichomonas_vaginalis")
    L(f"\n  T. vaginalis: {tvag_binds} single-primer 3'-binds, {tvag_amps} amplicons")
    L(f"  FULL off-target panel: {len(amp_rows)} amplicons")
    return len(amp_rows), tvag_amps


if __name__ == "__main__":
    n_all, n_tvag = main()
    logp = PROJ / "v5" / "00_logs" / "t5_offtarget.log"
    logp.write_text("\n".join(log) + "\n")
    print(f"Log -> {logp}")
