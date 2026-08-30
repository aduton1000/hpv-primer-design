#!/usr/bin/env python
"""
v4 TASK 6 - genome coverage incl. origin per type, recomputed on the v4 design
(HPV58 amp3 moved; confirm HPV58 still 100% incl origin; HPV52 must stay 100%).
Coverage = union of amplicon spans (left primer start -> right primer end), circular,
over the canonical genome length; origin = canonical position 0 covered.
"""
import json
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
DES = PROJ / "v5" / "04_design"
OUT = PROJ / "v5" / "06_reports"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
log = []
def L(m): print(m); log.append(str(m))


def coverage_circular(intervals, Lc):
    pts = [False] * Lc
    for a, b in intervals:
        a %= Lc; b %= Lc
        if a <= b:
            for x in range(a, b): pts[x] = True
        else:
            for x in range(a, Lc): pts[x] = True
            for x in range(0, b): pts[x] = True
    return round(100.0 * sum(pts) / Lc, 2), pts[0]


def main():
    L("=== T6 coverage incl. origin (v4) ===")
    rows = ["type\tgenome_len\tn_amplicons\tcoverage_pct\torigin_covered\tn_junction\tmin_overlap_bp"]
    allcov = True
    for t in TYPES:
        sch = json.load(open(DES / f"HPV{t}.scheme.json"))
        Lc = int(sch["genome_len"])
        ivs = []
        njunc = 0
        spans = []
        for ln in (DES / f"HPV{t}.amplicons.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t")
            s, e = int(c[3]), int(c[4])
            ivs.append((s, e)); spans.append((s, e))
            if c[5] == "Y": njunc += 1
        cov, origin = coverage_circular(ivs, Lc)
        # min overlap among overlapping amplicon pairs
        def sset(s, e):
            return set(range(s, e)) if s <= e else set(range(s, Lc)) | set(range(0, e))
        ss = [sset(s, e) for s, e in spans]
        ovs = [len(ss[i] & ss[j]) for i in range(len(ss)) for j in range(i+1, len(ss)) if ss[i] & ss[j]]
        min_ov = min(ovs) if ovs else 0
        rows.append(f"HPV{t}\t{Lc}\t{len(ivs)}\t{cov}\t{'Y' if origin else 'N'}\t{njunc}\t{min_ov}")
        L(f"HPV{t}: len={Lc} amps={len(ivs)} coverage={cov}% origin={'Y' if origin else 'N'} "
          f"junctions={njunc} min_overlap={min_ov}bp")
        if cov < 100.0 or not origin:
            allcov = False
    (OUT / "coverage_map.tsv").write_text("\n".join(rows) + "\n")
    L(f"\nAll types 100% coverage incl origin: {allcov}")
    return allcov


if __name__ == "__main__":
    main()
    (PROJ / "v5" / "00_logs" / "t5_coverage.log").write_text("\n".join(log) + "\n")
    print(f"Log -> {PROJ/'v4/00_logs/t6_coverage.log'}")
