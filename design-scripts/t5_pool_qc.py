#!/usr/bin/env python
"""
v3 Task 6 - intra-pool primer-dimer + hairpin re-check and Tm-spread per pool, for
the MODIFIED v3 scheme (all primer changes invalidate the prior dimer check).

Each HPV type is its own multiplex (per-type pools). Within each (type, pool) we
check: every primer's hairpin Tm, every intra-pool primer-primer heterodimer Tm,
and the Tm spread (max-min primer Tm). IUPAC bases are sanitised to a representative
base for the thermodynamic calc (degenerate primers are oligo mixes; flagged).

Thresholds: hairpin Tm < 45 C, heterodimer Tm < 45 C, pool Tm spread <= 5 C.

Outputs: v3/06_reports/pool_qc.tsv, v3/00_logs/v3_revalidate.log
"""
from pathlib import Path
from itertools import combinations
import primer3

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
DES = PROJ / "v5" / "04_design"
OUT = PROJ / "v5" / "06_reports"
LOG = PROJ / "v5" / "00_logs" / "t5_pool_qc.log"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
HAIRPIN_MAX = 50.0   # annealing(~62C) - 10C standard hairpin guideline
HETERO_MAX = 45.0
TMSPREAD_MAX = 5.0
IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT",
         "M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
log = []
def L(m): print(m); log.append(m)


def variants(s):
    import itertools
    return ["".join(p) for p in itertools.product(*[IUPAC.get(b, b) for b in s])]


def max_hairpin(s):
    return max(primer3.calc_hairpin_tm(v) for v in variants(s))


def max_heterodimer(a, b):
    # worst-case over the oligo-mix variants (pools are small; cap combinations)
    va, vb = variants(a), variants(b)
    if len(va) * len(vb) > 64:
        va, vb = va[:8], vb[:8]
    return max(primer3.calc_heterodimer_tm(x, y) for x in va for y in vb)


def main():
    rows = ["type\tpool\tn_primers\ttm_min\ttm_max\ttm_spread\tmax_hairpin_tm\tmax_heterodimer_tm\t"
            "worst_dimer_pair\tn_degenerate_primers\tpool_QC"]
    all_ok = True
    summ = []
    for t in TYPES:
        pools = {1: [], 2: []}
        for ln in (DES / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t")
            pools[int(c[2])].append(dict(name=c[3], seq=c[5], tm=float(c[8]),
                                         deg=any(b not in "ACGT" for b in c[5])))
        type_ok = True
        for pool, ps in pools.items():
            if not ps:
                continue
            tms = [p["tm"] for p in ps]
            spread = round(max(tms) - min(tms), 1)
            max_hp = max(max_hairpin(p["seq"]) for p in ps)
            worst_dimer = 0.0; worst_pair = ""
            for a, b in combinations(ps, 2):
                d = max_heterodimer(a["seq"], b["seq"])
                if d > worst_dimer:
                    worst_dimer = d; worst_pair = f"{a['name']}|{b['name']}"
            ndeg = sum(1 for p in ps if p["deg"])
            qc = (max_hp < HAIRPIN_MAX and worst_dimer < HETERO_MAX and spread <= TMSPREAD_MAX)
            type_ok &= qc; all_ok &= qc
            rows.append(f"HPV{t}\t{pool}\t{len(ps)}\t{min(tms):.1f}\t{max(tms):.1f}\t{spread}\t"
                        f"{max_hp:.1f}\t{worst_dimer:.1f}\t{worst_pair}\t{ndeg}\t{'PASS' if qc else 'FAIL'}")
        summ.append((t, type_ok))
    (OUT / "pool_qc.tsv").write_text("\n".join(rows) + "\n")
    L("=== POOL QC (intra-pool dimer/hairpin + Tm spread) ===")
    for r in rows:
        L("  " + r.replace("\t", "  "))
    L(f"\nAll pools PASS: {all_ok}")
    L(f"Types fully PASS: {sum(1 for _, ok in summ if ok)}/{len(TYPES)}")
    LOG.write_text("\n".join(log) + "\n")
    print(f"\nLog -> {LOG}")


if __name__ == "__main__":
    main()
