#!/usr/bin/env python
"""
v2 STAGE 6 - inclusivity validation (in-silico PCR of the final scheme against the
FULL per-type panel), broken down BY LINEAGE, with held-out stratification.

In-silico PCR (IUPAC-aware, circular via doubled genome):
  * a primer PRIMES a genome if its 3' terminal 5 nt match EXACTLY and the whole
    primer matches with <=3 total mismatches;
  * an amplicon AMPLIFIES a genome if its F primer (plus strand) and R primer
    (revcomp on plus strand) prime in convergent orientation with product 2-4 kb.
Reported per amplicon: % of panel genomes amplified, split by lineage and by
African-country origin. Strata: DESIGN set (deduped canonical, used for design),
HELD-OUT (near-duplicates removed in dedup -> unseen), DIVERGENT (set-aside).

Outputs:
  v2/05_specificity/inclusivity_by_amplicon.tsv
  v2/05_specificity/inclusivity_by_lineage.tsv
  v2/05_specificity/inclusivity_by_stratum.tsv
  v2/06_reports/inclusivity_summary.tsv
  v2/00_logs/v2_08_inclusivity.log
"""
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
RAW = PROJ / "v2" / "02_data" / "raw"
DD = PROJ / "v2" / "02_data" / "dedup"
DES = PROJ / "v5" / "04_design"
SPEC = PROJ / "v5" / "05_specificity"
LOG = PROJ / "v5" / "00_logs" / "t5_inclusivity.log"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
AMIN, AMAX = 2000, 4000          # design target
# In-silico "does it amplify?" uses a PCR-realistic window, not the design target:
# indels shift the actual product a little, and any convergent product a polymerase
# can make counts as amplification (a 1961 bp product still amplifies fine).
PCR_MIN, PCR_MAX = 600, 4800   # 600 accommodates short single-type gap amplicons (HPV52 L1 ~825-1117bp)
MAXMISM = 3
IUPAC = {"A":set("A"),"C":set("C"),"G":set("G"),"T":set("T"),"R":set("AG"),"Y":set("CT"),
         "S":set("GC"),"W":set("AT"),"K":set("GT"),"M":set("AC"),"B":set("CGT"),"D":set("AGT"),
         "H":set("ACT"),"V":set("ACG"),"N":set("ACGT")}
AFR = ["ghana","nigeria","south africa","zimbabwe","togo","kenya","uganda","tanzania","rwanda",
       "gambia","senegal","mozambique","zambia","cameroon","ethiopia","mali","burkina","ivoire",
       "benin","botswana","malawi","namibia","angola","congo","sudan","madagascar","guinea"]
log = []
def L(m): print(m); log.append(m)


def match_at(pattern, seq, pos, anchor):
    """Does pattern match seq starting at pos (IUPAC on pattern), with 3' anchor
    (anchor='end' -> last5 exact; 'start' -> first5 exact) and <=MAXMISM total?"""
    if pos < 0 or pos + len(pattern) > len(seq):
        return False
    mism = 0
    for k, pb in enumerate(pattern):
        sb = seq[pos + k]
        if sb not in IUPAC.get(pb, set()):
            mism += 1
            if mism > MAXMISM:
                return False
    # 3' clamp exact (no mismatch in the anchored 5)
    rng = range(len(pattern) - 5, len(pattern)) if anchor == "end" else range(0, 5)
    for k in rng:
        if seq[pos + k] not in IUPAC.get(pattern[k], set()):
            return False
    return True


def find_sites(pattern, seq, anchor):
    """All start positions where pattern primes seq (exact 5-mer seed + verify)."""
    seed = pattern[-5:] if anchor == "end" else pattern[:5]
    # seed has no IUPAC (3' clamp is mono by design) -> exact string search
    offset = len(pattern) - 5 if anchor == "end" else 0
    sites, start = [], 0
    while True:
        i = seq.find(seed, start)
        if i < 0: break
        pos = i - offset
        if match_at(pattern, seq, pos, anchor):
            sites.append(pos)
        start = i + 1
    return sites


def _amp_one_strand(fseq, rseq, dbl, glen):
    f_sites = find_sites(fseq, dbl, "end")                        # F on plus strand
    rpat = str(Seq(rseq).reverse_complement())                   # R binding region on plus strand
    r_sites = find_sites(rpat, dbl, "start")                     # R 3' anchors at its start
    for pf in f_sites:
        for pr in r_sites:
            prod = (pr + len(rpat)) - pf
            if PCR_MIN <= prod <= PCR_MAX and pf < pr + len(rpat) and pf < glen:
                return True
    return False


def amplifies(fseq, rseq, genome_doubled, glen, genome_rc_doubled):
    """True if F+R make a 2-4 kb convergent product on EITHER strand of the genome
    (deposited genomes vary in strand/rotation vs the canonical primer frame)."""
    return (_amp_one_strand(fseq, rseq, genome_doubled, glen)
            or _amp_one_strand(fseq, rseq, genome_rc_doubled, glen))


def main():
    # lineage + country per accession
    lineage, country = {}, {}
    for ln in (PROJ / "v2" / "02_data" / "lineage_assignments.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); lineage[c[0]] = c[2]; country[c[0]] = c[5]
    for t in TYPES:                       # also country for genomes not in lineage file
        for ln in (RAW / f"HPV{t}.provenance.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t"); country.setdefault(c[0], c[4])
    def is_afr(acc): return any(a in (country.get(acc, "") or "").lower() for a in AFR)

    amp_rows, lin_rows, strat_rows = [], [], []
    summary = []
    for t in TYPES:
        canon_accs = {r.id.split("|")[1] for r in SeqIO.parse(str(DD / f"HPV{t}.norm.fasta"), "fasta")}
        test_accs = set((DD / f"HPV{t}.test_accessions.txt").read_text().split()) \
            if (DD / f"HPV{t}.test_accessions.txt").exists() else set()
        design_accs = canon_accs - test_accs          # 80% train = design conservation set
        div_accs = {r.id.split("|")[1] for r in SeqIO.parse(str(DD / f"HPV{t}.divergent.fasta"), "fasta")} \
            if (DD / f"HPV{t}.divergent.fasta").exists() else set()
        # amplicons (F/R seqs) for this type
        amps = defaultdict(dict)
        for ln in (DES / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t"); amps[c[1]][c[4]] = c[5]
        panel = list(SeqIO.parse(str(RAW / f"HPV{t}.gb.fasta"), "fasta"))
        # precompute doubled genomes
        genomes = []
        for r in panel:
            s = str(r.seq).upper(); rcs = str(Seq(s).reverse_complement())
            genomes.append((r.id.split("|")[1], s + s, len(s), rcs + rcs))

        per_amp_amp = {}        # amp -> set(acc amplified)
        for amp, fr in amps.items():
            fseq, rseq = fr["F"], fr["R"]
            amped = set()
            for acc, dbl, glen, rcdbl in genomes:
                if amplifies(fseq, rseq, dbl, glen, rcdbl):
                    amped.add(acc)
            per_amp_amp[amp] = amped

        # genome amplified by ALL amplicons?
        all_accs = [g[0] for g in genomes]
        full = {acc for acc in all_accs if all(acc in per_amp_amp[a] for a in amps)}

        # strata
        def stratum(acc):
            if acc in design_accs: return "design"
            if acc in div_accs: return "divergent"
            return "heldout"
        # per amplicon inclusivity (overall + african)
        nA = sum(1 for acc in all_accs if is_afr(acc))
        for amp in sorted(amps, key=int):
            amped = per_amp_amp[amp]
            inc = 100.0 * len(amped) / len(all_accs)
            afr_amped = sum(1 for acc in amped if is_afr(acc))
            inc_afr = 100.0 * afr_amped / nA if nA else 100.0
            amp_rows.append((f"HPV{t}", amp, len(all_accs), round(inc, 1), nA, round(inc_afr, 1),
                             "FLAG" if (inc < 95 or (nA and inc_afr < 95)) else "OK"))
        # per lineage (using 'full' = amplified by all amplicons = genome fully tiled)
        by_lin = defaultdict(lambda: [0, 0])
        for acc in all_accs:
            lin = lineage.get(acc, "NA")
            by_lin[lin][1] += 1
            if acc in full: by_lin[lin][0] += 1
        for lin, (a, n) in sorted(by_lin.items()):
            lin_rows.append((f"HPV{t}", lin, n, a, round(100.0 * a / n, 1) if n else 0))
        # per stratum
        by_str = defaultdict(lambda: [0, 0])
        for acc in all_accs:
            st = stratum(acc); by_str[st][1] += 1
            if acc in full: by_str[st][0] += 1
        for st in ("design", "heldout", "divergent"):
            a, n = by_str[st]
            if n: strat_rows.append((f"HPV{t}", st, n, a, round(100.0 * a / n, 1)))
        worst_amp = min((100.0 * len(per_amp_amp[a]) / len(all_accs)) for a in amps)
        min_afr = min(((100.0 * sum(1 for acc in per_amp_amp[a] if is_afr(acc)) / nA) if nA else 100)
                      for a in amps)
        summary.append((t, len(all_accs), nA, round(100.0 * len(full) / len(all_accs), 1),
                        round(worst_amp, 1), round(min_afr, 1)))
        L(f"HPV{t}: panel={len(all_accs)} (African {nA}); fully-tiled={100.0*len(full)/len(all_accs):.1f}% "
          f"worst-amplicon={worst_amp:.1f}% min-African-amplicon={min_afr:.1f}%")

    with open(SPEC / "inclusivity_by_amplicon.tsv", "w") as fh:
        fh.write("type\tamplicon\tn_panel\tinclusivity_pct\tn_african\tafrican_inclusivity_pct\tflag\n")
        for r in amp_rows: fh.write("\t".join(str(x) for x in r) + "\n")
    with open(SPEC / "inclusivity_by_lineage.tsv", "w") as fh:
        fh.write("type\tlineage\tn\tfully_amplified\tpct\n")
        for r in lin_rows: fh.write("\t".join(str(x) for x in r) + "\n")
    with open(SPEC / "inclusivity_by_stratum.tsv", "w") as fh:
        fh.write("type\tstratum\tn\tfully_amplified\tpct\n")
        for r in strat_rows: fh.write("\t".join(str(x) for x in r) + "\n")
    with open(PROJ / "v5" / "06_reports" / "inclusivity_summary.tsv", "w") as fh:
        fh.write("type\tn_panel\tn_african\tfully_tiled_pct\tworst_amplicon_pct\tmin_african_amplicon_pct\n")
        for t, n, na, ft, wa, ma in summary:
            fh.write(f"HPV{t}\t{n}\t{na}\t{ft}\t{wa}\t{ma}\n")

    L("\n=== INCLUSIVITY SUMMARY ===")
    L(f"{'type':7}{'panel':>6}{'afr':>5}{'fulltiled%':>11}{'worstAmp%':>10}{'minAfrAmp%':>11}")
    for t, n, na, ft, wa, ma in summary:
        L(f"HPV{t:<4}{n:>6}{na:>5}{ft:>11}{wa:>10}{ma:>11}")
    LOG.write_text("\n".join(log) + "\n")
    print(f"\nLog -> {LOG}")


if __name__ == "__main__":
    main()
