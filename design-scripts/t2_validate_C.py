#!/usr/bin/env python
"""
v5 TASK 2 - in-silico PCR all 6 HPV58 amplicons (v4 primers, incl redesigned amp3)
against the complete lineage-C genomes from T1. Per amplicon per genome: both-primer
3'-anchored binding, 3'-clamp exact match on C, total mismatch count + positions,
predicted product, PASS/FAIL. Summarised per amplicon over all C genomes.

Output: v5/05_specificity/hpv58_C_validation.tsv
"""
import importlib.util
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
RAW = PROJ / "v2" / "02_data" / "raw"
DESV4 = PROJ / "v4" / "04_design"
PROVC = PROJ / "v5" / "02_data" / "hpv58_C_provenance.tsv"
OUT = PROJ / "v5" / "05_specificity"; OUT.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("v3i", str(PROJ / "v3/scripts/v3_inclusivity.py"))
iz = importlib.util.module_from_spec(spec); spec.loader.exec_module(iz)
IUP = iz.IUPAC
PCR_MIN, PCR_MAX = 600, 4800
log = []
def L(m): print(m); log.append(str(m))


def prime_best(primer, seq, anchor):
    """Return (primes_bool, best_total_mm, clamp5_mm, mismatch_positions_3p_indexed).
    3'-anchored prime requires clamp exact (5nt) AND total mm<=3."""
    L_ = len(primer); best = None
    for strand, s in (("+", seq + seq), ("-", str(Seq(seq).reverse_complement()) * 2)):
        for i in range(len(s) - L_):
            seg = s[i:i + L_]
            mm = sum(1 for a, b in zip(primer, seg) if b not in IUP.get(a, {a}))
            if best is None or mm < best[0]:
                rng = range(L_ - 5, L_) if anchor == "end" else range(0, 5)
                clampmm = sum(1 for k in rng if seg[k] not in IUP.get(primer[k], {primer[k]}))
                # mismatch positions counted from the 3' end (1 = 3'-terminal base)
                mmpos = []
                for k in range(L_):
                    if seg[k] not in IUP.get(primer[k], {primer[k]}):
                        mmpos.append(L_ - k if anchor == "end" else k + 1)
                best = (mm, clampmm, sorted(mmpos), strand, i)
    mm, clampmm, mmpos, strand, ipos = best
    primes = (mm <= 3 and clampmm == 0)
    return primes, mm, clampmm, mmpos, strand, ipos


def main():
    amps = {}
    for ln in (DESV4 / "HPV58.primers.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); amps.setdefault(c[1], {})[c[4]] = c[5]
    products = {}
    for ln in (DESV4 / "HPV58.amplicons.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); products[c[0]] = int(c[2])
    seqs = {r.id.split("|")[1]: str(r.seq).upper()
            for r in SeqIO.parse(str(RAW / "HPV58.gb.fasta"), "fasta")}
    Cgen = [ln.split("\t")[0] for ln in PROVC.read_text().splitlines()[1:]
            if ln.split("\t")[-1] == "YES"]
    L(f"=== T2: HPV58 amplicon validation vs {len(Cgen)} complete lineage-C genomes ===")
    L(f"C genomes: {Cgen}")

    rows = ["amplicon\tn_C_genomes\tn_PASS\tn_FAIL\tworst_left_mm\tworst_right_mm\t"
            "left_clamp_ok_all\tright_clamp_ok_all\tproduct_range_bp\tverdict\tfailing_C_accessions"]
    per_amp_pass = {}
    for amp in ("1", "2", "3", "4", "5", "6"):
        F, R = amps[amp]["F"], amps[amp]["R"]
        rpat = str(Seq(R).reverse_complement())
        npass = 0; fails = []
        wlmm = 0; wrmm = 0; lclamp_all = True; rclamp_all = True; prods = []
        for acc in Cgen:
            s = seqs[acc]
            rcs = str(Seq(s).reverse_complement())
            amp_ok = iz.amplifies(F, R, s + s, len(s), rcs + rcs)
            lp, lmm, lcl, lpos, lstr, lipos = prime_best(F, s, "end")
            rp, rmm, rcl, rpos, rstr, ripos = prime_best(rpat, s, "start")
            wlmm = max(wlmm, lmm); wrmm = max(wrmm, rmm)
            lclamp_all &= (lcl == 0); rclamp_all &= (rcl == 0)
            if amp_ok:
                npass += 1
                # product size
                for tmpl in (s + s, rcs + rcs):
                    fs = iz.find_sites(F, tmpl, "end"); rs = iz.find_sites(rpat, tmpl, "start")
                    for pf in fs:
                        for pr in rs:
                            prod = (pr + len(rpat)) - pf
                            if PCR_MIN <= prod <= PCR_MAX and pf < pr + len(rpat) and pf < len(s):
                                prods.append(prod)
            else:
                fails.append(acc)
        per_amp_pass[amp] = npass
        verdict = "PASS" if npass == len(Cgen) else ("PARTIAL" if npass else "FAIL")
        prng = f"{min(prods)}-{max(prods)}" if prods else "NA"
        rows.append(f"{amp}\t{len(Cgen)}\t{npass}\t{len(fails)}\t{wlmm}\t{wrmm}\t"
                    f"{'Y' if lclamp_all else 'N'}\t{'Y' if rclamp_all else 'N'}\t{prng}\t"
                    f"{verdict}\t{','.join(fails) if fails else '-'}")
        L(f"  amp{amp}: {npass}/{len(Cgen)} PASS | worst LEFT mm={wlmm} RIGHT mm={wrmm} | "
          f"clamps L:{'ok' if lclamp_all else 'X'} R:{'ok' if rclamp_all else 'X'} | "
          f"product {prng} | {verdict}{'  FAILS: '+','.join(fails) if fails else ''}")
    (OUT / "hpv58_C_validation.tsv").write_text("\n".join(rows) + "\n")
    allpass = all(per_amp_pass[a] == len(Cgen) for a in per_amp_pass)
    L(f"\n=== VERDICT: all 6 HPV58 amplicons amplify ALL {len(Cgen)} complete C genomes: "
      f"{allpass} ===")
    L(f"wrote {OUT/'hpv58_C_validation.tsv'}")
    return allpass, per_amp_pass, len(Cgen)


if __name__ == "__main__":
    allpass, per, n = main()
    (PROJ / "v5" / "00_logs" / "t2_validate_C.log").write_text("\n".join(log) + "\n")
    print(f"Log -> {PROJ/'v5/00_logs/t2_validate_C.log'}")
