#!/usr/bin/env python
"""
v5 TASK 4 - light audit for the HPV58-C trap in the other 11 types.

For each type, group panel genomes by major lineage (A/B/C/D from lineage_hint;
'' -> NA). For each NON-A lineage present, report:
  n_total, n_complete (non-N at ALL amp binding sites), n_complete_passing_all_amps,
  and a status:
    VALIDATED      - >=1 complete genome of that lineage is amplified by every amplicon
    BLIND_SPOT_Ngap- lineage present but EVERY representative is N-gapped at a binding
                     site (the HPV58-C v4 trap: present but never truly tested)
    INCLUSIVITY_GAP- complete representative(s) exist but at least one fails an amplicon
Annotates African-relevant non-A lineages per HPV lineage literature.

This is verification/reporting: flag blind spots; do not redesign unless a concrete
untested-lineage gap is found.

Output: v5/05_specificity/lineage_audit.tsv
"""
import importlib.util
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
RAW = PROJ / "v2" / "02_data" / "raw"
DES = PROJ / "v5" / "04_design"
OUT = PROJ / "v5" / "05_specificity"
spec = importlib.util.spec_from_file_location("v3i", str(PROJ / "v3/scripts/v3_inclusivity.py"))
iz = importlib.util.module_from_spec(spec); spec.loader.exec_module(iz)
IUP = iz.IUPAC
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 59, 66]   # the OTHER 11 (HPV58 done in T1/T2)
PCR_MIN, PCR_MAX = 600, 4800
# African-relevant non-A lineages per type (HPV lineage literature; B/C are the
# classic "African" branches for the alpha-9/7 high-risk types).
AFRICAN_NONA = {
    16: {"B", "C"},   # B=African-1, C=African-2 (D is AA/As-Am, non-African)
    18: {"B", "C"},   # B/C African branches
    31: {"B", "C"},
    33: {"B", "C"},
    35: {"B"},        # HPV35 mostly A; B non-A
    39: {"B", "C"},
    45: {"B", "C"},
    52: {"B", "C", "D"},
    56: {"B"},
    59: {"B"},
    66: {"B"},
}
log = []
def L(m): print(m); log.append(str(m))


def prime_pos(primer, seq, anchor):
    out = []; L_ = len(primer)
    for strand, s in (("+", seq + seq), ("-", str(Seq(seq).reverse_complement()) * 2)):
        for i in range(len(s) - L_):
            seg = s[i:i + L_]
            rng = range(L_ - 5, L_) if anchor == "end" else range(0, 5)
            if all(seg[k] in IUP.get(primer[k], {primer[k]}) for k in rng):
                mm = sum(1 for a, b in zip(primer, seg) if b not in IUP.get(a, {a}))
                if mm <= 3:
                    out.append((i, strand))
    return out


def site_complete(F, R, seq, product):
    rpat = str(Seq(R).reverse_complement())
    fp = prime_pos(F, seq, "end"); rp = prime_pos(rpat, seq, "start")
    def winN(strand, start, length):
        s = (seq + seq) if strand == "+" else str(Seq(seq).reverse_complement()) * 2
        start = max(0, start); w = s[start:start + length]
        return sum(1 for c in w if c not in "ACGT"), len(w)
    if fp and rp:
        return True
    if fp and not rp:
        i, strand = fp[0]; f3 = i + len(F)
        exp = f3 + product - len(F) - len(rpat); nN, wl = winN(strand, exp - 3, len(rpat) + 6)
        return nN < 0.5 * wl
    if rp and not fp:
        i, strand = rp[0]; exp = i + len(rpat) - product + len(F)
        nN, wl = winN(strand, exp - len(F) - 3, len(F) + 6)
        return nN < 0.5 * wl
    return sum(1 for c in seq if c not in "ACGT") < 0.05 * len(seq)


def main():
    L("=== T4: African non-A lineage audit (HPV58-C trap check) for the other 11 types ===")
    rows = ["type\tlineage\tafrican_relevant\tn_total\tn_complete_at_sites\tn_complete_pass_all_amps\t"
            "status"]
    flagged = []
    for t in TYPES:
        amps = {}
        for ln in (DES / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t"); amps.setdefault(c[1], {})[c[4]] = c[5]
        products = {}
        for ln in (DES / f"HPV{t}.amplicons.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t"); products[c[0]] = int(c[2])
        prov = {c[0]: c for c in (ln.split("\t") for ln in
                (RAW / f"HPV{t}.provenance.tsv").read_text().splitlines()[1:])}
        seqs = {r.id.split("|")[1]: str(r.seq).upper()
                for r in SeqIO.parse(str(RAW / f"HPV{t}.gb.fasta"), "fasta")}
        # group by major lineage
        bylin = defaultdict(list)
        for acc in seqs:
            lh = prov.get(acc, [""]*7)[6] if len(prov.get(acc, [])) > 6 else ""
            major = lh[0] if lh else "NA"
            bylin[major].append(acc)
        nonA = sorted([k for k in bylin if k not in ("A", "NA")])
        for lin in nonA:
            accs = bylin[lin]
            ncomplete = 0; npass = 0
            for acc in accs:
                s = seqs[acc]
                comp = all(site_complete(amps[a]["F"], amps[a]["R"], s, products[a])
                           for a in ("1", "2", "3", "4", "5", "6") if a in amps) \
                       if all(x in amps for x in ("1","2","3","4","5","6")) else \
                       all(site_complete(amps[a]["F"], amps[a]["R"], s, products[a]) for a in amps)
                if comp:
                    ncomplete += 1
                    rcs = str(Seq(s).reverse_complement())
                    allamp = all(iz.amplifies(amps[a]["F"], amps[a]["R"], s + s, len(s), rcs + rcs)
                                 for a in amps)
                    if allamp:
                        npass += 1
            afr = lin in AFRICAN_NONA.get(t, set())
            if ncomplete == 0:
                status = "BLIND_SPOT_Ngap"
            elif npass >= 1:
                status = "VALIDATED"
            else:
                status = "INCLUSIVITY_GAP"
            if afr and status != "VALIDATED":
                flagged.append((f"HPV{t}", lin, status, len(accs), ncomplete))
            rows.append(f"HPV{t}\t{lin}\t{'Y' if afr else 'n'}\t{len(accs)}\t{ncomplete}\t{npass}\t{status}")
            L(f"  HPV{t} lineage {lin} (African={'Y' if afr else 'n'}): total={len(accs)} "
              f"complete={ncomplete} pass_all={npass} -> {status}")
        if not nonA:
            rows.append(f"HPV{t}\t(none non-A)\tn\t0\t0\t0\tNO_NONA_LINEAGE")
            L(f"  HPV{t}: no non-A lineage in panel")
    (OUT / "lineage_audit.tsv").write_text("\n".join(rows) + "\n")
    L("\n=== BLIND-SPOT / GAP flags (African-relevant non-A not VALIDATED) ===")
    if flagged:
        for f in flagged:
            L(f"  FLAG: {f[0]} lineage {f[1]} -> {f[2]} (total {f[3]}, complete {f[4]})")
    else:
        L("  NONE - every African-relevant non-A lineage present has a complete,")
        L("  amplicon-validated representative. No HPV58-C-style trap in the other types.")
    L(f"wrote {OUT/'lineage_audit.tsv'}")
    return flagged


if __name__ == "__main__":
    flagged = main()
    (PROJ / "v5" / "00_logs" / "t4_lineage_audit.log").write_text("\n".join(log) + "\n")
    print(f"Log -> {PROJ/'v5/00_logs/t4_lineage_audit.log'}")
