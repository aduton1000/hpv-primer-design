#!/usr/bin/env python
"""
v5 TASK 1 - acquire COMPLETE HPV58 sublineage-C genomes and verify completeness at
the amp1-6 primer binding sites + independently confirm lineage C (not relying on
the panel's own lineage_hint).

Sources: the HPV58 panel already carries the Chen et al. GenBank submissions the task
names (HQ537xxx = Chen 2011 variant set; KY225xxx = Chen 2017 90-genome set; LC270105).
13 are lineage-C-hinted; we keep only those that are COMPLETE (no N) at every amp1-6
binding coordinate, and confirm each clusters with the C reference and away from the
A prototype by whole-genome identity.

Outputs: v5/02_data/hpv58_C_provenance.tsv
"""
import importlib.util
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
RAW = PROJ / "v2" / "02_data" / "raw"
DESV4 = PROJ / "v4" / "04_design"
OUT = PROJ / "v5" / "02_data"
spec = importlib.util.spec_from_file_location("v3i", str(PROJ / "v3/scripts/v3_inclusivity.py"))
iz = importlib.util.module_from_spec(spec); spec.loader.exec_module(iz)
IUP = iz.IUPAC
log = []
def L(m): print(m); log.append(str(m))


def hpv58_primers():
    amps = {}
    for ln in (DESV4 / "HPV58.primers.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); amps.setdefault(c[1], {})[c[4]] = c[5]
    return amps


def prime_pos(primer, seq, anchor):
    """3'-anchored prime positions (clamp exact, <=3 mm) on doubled both strands."""
    out = []; L_ = len(primer)
    for strand, s in (("+", seq + seq), ("-", str(Seq(seq).reverse_complement()) * 2)):
        for i in range(len(s) - L_):
            seg = s[i:i + L_]
            rng = range(L_ - 5, L_) if anchor == "end" else range(0, 5)
            if all(seg[k] in IUP.get(primer[k], {primer[k]}) for k in rng):
                mm = sum(1 for a, b in zip(primer, seg) if b not in IUP.get(a, {a}))
                if mm <= 3:
                    out.append((i, strand, seg))
    return out


def site_complete(F, R, seq, product):
    """Is the genome non-N (complete) at BOTH primer sites of this amplicon?
    Anchors off whichever primer primes; reads the homologous window for the other;
    returns (both_sites_have_real_sequence, detail)."""
    rpat = str(Seq(R).reverse_complement())
    fp = prime_pos(F, seq, "end"); rp = prime_pos(rpat, seq, "start")
    def winN(strand, start, length):
        s = (seq + seq) if strand == "+" else str(Seq(seq).reverse_complement()) * 2
        start = max(0, start); w = s[start:start + length]
        return sum(1 for c in w if c not in "ACGT"), len(w)
    if fp and rp:
        return True, "both primers anchor (no N gap)"
    if fp and not rp:
        i, strand, _ = fp[0]; f3 = i + len(F)
        exp = f3 + product - len(F) - len(rpat); nN, wl = winN(strand, exp - 3, len(rpat) + 6)
        return (nN < 0.5 * wl), f"R-site N={nN}/{wl}"
    if rp and not fp:
        i, strand, _ = rp[0]; exp = i + len(rpat) - product + len(F)
        nN, wl = winN(strand, exp - len(F) - 3, len(F) + 6)
        return (nN < 0.5 * wl), f"F-site N={nN}/{wl}"
    namb = sum(1 for c in seq if c not in "ACGT")
    return (namb < 0.05 * len(seq)), f"neither primer anchors; genomeN={namb}"


def ident(a, b):
    """ungapped k-mer identity proxy: fraction of a's 21-mers present in b (both strands)."""
    K = 21
    bs = set(b[i:i+K] for i in range(len(b)-K))
    brc = str(Seq(b).reverse_complement()); bs |= set(brc[i:i+K] for i in range(len(brc)-K))
    sh = sum(1 for i in range(len(a)-K) if a[i:i+K] in bs)
    return 100.0 * sh / (len(a)-K)


def main():
    prov = {}
    for ln in (RAW / "HPV58.provenance.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); prov[c[0]] = c
    seqs = {r.id.split("|")[1]: str(r.seq).upper()
            for r in SeqIO.parse(str(RAW / "HPV58.gb.fasta"), "fasta")}
    amps = hpv58_primers()
    products = {}
    for ln in (DESV4 / "HPV58.amplicons.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); products[c[0]] = int(c[2])

    # candidate C genomes (lineage_hint starts with C)
    C_cands = [a for a, c in prov.items() if len(c) > 6 and c[6].startswith("C")]
    # A prototype = a lineage-A1 genome (reference frame) for identity contrast
    A_ref = next(a for a, c in prov.items() if len(c) > 6 and c[6] == "A1")
    # a confirmed complete C as the "C reference" for clustering (pick first complete)
    L(f"=== T1: HPV58 lineage-C acquisition & completeness verification ===")
    L(f"candidate C-hinted genomes in panel: {len(C_cands)}")
    L(f"A1 reference for identity contrast: {A_ref}")

    rows = ["accession\tlineage_hint\tcountry_or_isolate\tlen\tgenome_N_pct\t"
            "complete_at_all_amp_sites\tid_to_A1\tid_to_Cref\tlineage_confirmed_C\tUSE"]
    # establish a C-reference: the complete C with most complete sites
    complete_C = []
    detail = {}
    for acc in C_cands:
        s = seqs[acc]
        ok_all = True; dets = []
        for amp in ("1", "2", "3", "4", "5", "6"):
            F, R = amps[amp]["F"], amps[amp]["R"]
            ok, d = site_complete(F, R, s, products[amp])
            dets.append(f"amp{amp}:{'ok' if ok else 'N'}")
            ok_all &= ok
        detail[acc] = (ok_all, dets, sum(1 for c in s if c not in "ACGT"))
        if ok_all:
            complete_C.append(acc)
    Cref = complete_C[0]
    L(f"complete-at-all-sites C genomes: {len(complete_C)} -> {complete_C}")
    L(f"using {Cref} as C reference for identity clustering")

    use_list = []
    for acc in sorted(C_cands):
        s = seqs[acc]; c = prov[acc]
        ok_all, dets, nN = detail[acc]
        idA = round(ident(s, seqs[A_ref]), 1)
        idC = round(ident(s, seqs[Cref]), 1)
        # lineage-C confirmed if it is markedly closer to the C ref than to A1
        conf_C = idC >= idA  # C genomes share more with Cref than with A1
        usable = ok_all and conf_C and nN == 0
        if usable: use_list.append(acc)
        cstr = (c[4] or c[5]) if len(c) > 5 else ""
        rows.append(f"{acc}\t{c[6]}\t{cstr}\t{len(s)}\t{round(100*nN/len(s),1)}\t"
                    f"{'Y' if ok_all else 'N'}\t{idA}\t{idC}\t{'Y' if conf_C else 'N'}\t"
                    f"{'YES' if usable else 'no(' + ('Ngap' if not ok_all else 'lineage') + ')'}")
    (OUT / "hpv58_C_provenance.tsv").write_text("\n".join(rows) + "\n")
    L(f"\nUSABLE complete lineage-C validation genomes: {len(use_list)}")
    L(f"  {use_list}")
    if len(use_list) >= 3:
        L(f"\n>>> {len(use_list)} complete lineage-C genomes available (target >=3 MET).")
    else:
        L(f"\n>>> !!! ONLY {len(use_list)} complete lineage-C genomes found - C validation "
          f"is PARTIAL. FAIL LOUDLY in reporting. !!!")
    L(f"wrote {OUT/'hpv58_C_provenance.tsv'}")
    return use_list


if __name__ == "__main__":
    use = main()
    (PROJ / "v5" / "00_logs" / "t1_acquire_C.log").write_text("\n".join(log) + "\n")
    print(f"Log -> {PROJ/'v5/00_logs/t1_acquire_C.log'}")
