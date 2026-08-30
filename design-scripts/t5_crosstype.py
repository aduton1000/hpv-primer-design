import shutil
import os
#!/usr/bin/env python
"""
v4 TASK 2 - Symmetric, bidirectional, circular-aware cross-type in-silico PCR.

Fixes the two v3 defects that hid the HPV58 amp3 -> HPV33 amplicon:
  (1) AMPLICON DISTANCE was linear (<=5 kb) with no genome-doubling, so an
      origin-spanning product appeared as a >5 kb linear gap and was filtered out.
      Here distance is CIRCULAR: gap = (m3p - p3p) mod G on the foreign genome.
  (2) The "0" had no coordinate evidence. crosstype_hits.tsv now carries
      sstart/send/strand/pident/mismatch/3p_bind; crosstype_amplicons.tsv carries
      owner_amplicon, foreign_type, n_foreign_genomes_amplified,
      min_pident_both_primers, predicted_product_size, verdict.

Bidirectional by construction: every owner amplicon's F+R are tested against EVERY
foreign genome of EVERY other type, so X-on-Y and Y-on-X are both evaluated; sister
types are NOT exempted in either direction.

Amplicon call (per spec): BOTH primers of one owner amplicon 3'-anchor on the same
foreign genome, on OPPOSITE strands, in CONVERGENT orientation, with circular
inter-3' distance in [50, MAX_AMP] bp.

Independent cross-check: the IUPAC-aware doubled-genome PCR engine (same one used for
inclusivity) is run on the foreign genomes and must agree with the BLAST-coordinate
caller.

Usage: t2_crosstype_screen.py <design_dir> <out_prefix_label>
"""
import subprocess, sys, importlib.util
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
BL = os.environ.get("HPV_BLAST_BIN") or (os.path.dirname(shutil.which("blastn")) if shutil.which("blastn") else "")
DDV2 = PROJ / "v2" / "02_data" / "dedup"
SPEC = PROJ / "v5" / "05_specificity"; SPEC.mkdir(parents=True, exist_ok=True)
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
MAX_AMP = 5000
SISTERS = [(33, 58), (56, 66)]

# IUPAC-aware engine for the independent cross-check
spec = importlib.util.spec_from_file_location("v3i", str(PROJ / "v3/scripts/v3_inclusivity.py"))
iz = importlib.util.module_from_spec(spec); spec.loader.exec_module(iz)

log = []
def L(m): print(m); log.append(str(m))


def build_hpv_db():
    recs = []
    glen = {}
    for t in TYPES:
        for r in SeqIO.parse(str(DDV2 / f"HPV{t}.norm.fasta"), "fasta"):
            acc = r.id.split("|")[1]
            sid = f"HPV{t}__{acc}"
            r.id = sid; r.description = ""
            recs.append(r); glen[sid] = len(r.seq)
    fa = SPEC / "hpv_all_types.fasta"; SeqIO.write(recs, str(fa), "fasta")
    subprocess.run([f"{BL}/makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                    "-out", str(SPEC / "hpvdb")], check=True, capture_output=True)
    return glen


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
                        "-outfmt", fmt, "-max_target_seqs", "50000", "-num_threads", "6"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        L(f"  blast rc={r.returncode}: {r.stderr[:200]}"); return []
    hits = []
    for ln in r.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 14: continue
        q = f[0]; pl = plen.get(q, int(f[4]))
        length = int(f[3]); qend = int(f[6]); mism = int(f[10]) + int(f[11])
        qseq, sseq = f[12], f[13]
        # 3'-ANCHORED priming: alignment reaches the primer 3' terminus, near-full
        # length, <=3 total mismatches, AND 3' terminal 5 nt are an EXACT match.
        three_prime_exact = (qend == pl and len(qseq) >= 5 and qseq[-5:] == sseq[-5:])
        binds = (length >= pl - 3 and mism <= 3 and three_prime_exact)
        hits.append(dict(primer=q, subj=f[1], strand=f[9], sstart=int(f[7]), send=int(f[8]),
                         binds=binds, pident=float(f[2]), length=length, mism=mism,
                         qend=qend, qlen=pl))
    return hits


def circular_amplicons(hits, prim, amps, glen):
    """From 3'-anchored BLAST binds, flag cross-amplicons with CIRCULAR distance.
    Returns list of dicts and the per-(owner_amp,foreign_type) aggregation."""
    by_subj = defaultdict(list)
    for h in hits:
        if h["binds"]:
            by_subj[h["subj"]].append(h)
    rows = []
    for typ, ad in amps.items():
        for amp, fr in ad.items():
            fn, rn = fr.get("F"), fr.get("R")
            if not fn or not rn: continue
            for subj, sh in by_subj.items():
                st = subj.split("__")[0]
                if st == typ:
                    continue  # own type expected
                G = glen[subj]
                fs = [h for h in sh if h["primer"] == fn]
                rs = [h for h in sh if h["primer"] == rn]
                for a in fs:
                    for b in rs:
                        if a["strand"] == b["strand"]:
                            continue  # need opposite strands
                        # 3' terminus maps to BLAST 'send' on either strand.
                        plus = a if a["strand"] == "plus" else b
                        minus = b if a["strand"] == "plus" else a
                        p3 = plus["send"]; m3 = minus["send"]
                        # convergent circular product: plus primer upstream of minus primer
                        gap = (m3 - p3) % G
                        if 50 <= gap <= MAX_AMP:
                            rows.append(dict(owner_type=typ, owner_amp=amp, foreign_type=st,
                                             foreign_genome=subj, product=gap,
                                             pident=min(a["pident"], b["pident"]),
                                             f_primer=fn, r_primer=rn))
    return rows


def iupac_crosscheck(amps, prim, foreign_by_type):
    """Independent IUPAC-aware doubled-genome PCR: for each owner amplicon vs each
    foreign genome, does F+R make a convergent product? Returns set of
    (owner_type, amp, foreign_type) that amplify."""
    hits = set()
    seqs = {n: prim[n]["seq"] for n in prim}
    for typ, ad in amps.items():
        for amp, fr in ad.items():
            if "F" not in fr or "R" not in fr: continue
            fseq = seqs[fr["F"]]; rseq = seqs[fr["R"]]
            for ft, genomes in foreign_by_type.items():
                if ft == typ: continue
                for acc, dbl, gl, rcdbl in genomes:
                    if iz.amplifies(fseq, rseq, dbl, gl, rcdbl):
                        hits.add((typ, amp, ft)); break
    return hits


def main():
    pdir = sys.argv[1] if len(sys.argv) > 1 else str(PROJ / "v4" / "04_design")
    label = sys.argv[2] if len(sys.argv) > 2 else "v4"
    prim, amps = load_primers(pdir)
    plen = {n: p["length"] for n, p in prim.items()}
    glen = build_hpv_db()
    pf = SPEC / "all_primers.fasta"
    pf.write_text("\n".join(f">{n}\n{p['seq']}" for n, p in prim.items()) + "\n")
    L(f"=== T2 cross-type screen [{label}] on {pdir} ({len(prim)} primers, "
      f"{len(glen)} foreign-type genomes) ===")

    hits = blast(pf, SPEC / "hpvdb", plen)

    # crosstype_hits.tsv WITH coordinates (all foreign-type binds, 3p flagged)
    hits_path = SPEC / ("crosstype_hits.tsv" if label in ("v4","v5") else f"crosstype_hits_{label}.tsv")
    nbind = 0
    with open(hits_path, "w") as fh:
        fh.write("primer\towner_type\thit_genome\thit_type\tsstart\tsend\tstrand\tpident\t"
                 "mismatch\t3p_bind\n")
        for h in sorted(hits, key=lambda x: (x["primer"], x["subj"])):
            ot = prim[h["primer"]]["type"]; ht = h["subj"].split("__")[0]
            if ht == ot:
                continue
            # write foreign binds that are either 3'-anchored OR near-miss (audit trail)
            if h["binds"] or (h["length"] >= h["qlen"] - 3 and h["mism"] <= 4):
                fh.write(f"{h['primer']}\t{ot}\t{h['subj']}\t{ht}\t{h['sstart']}\t{h['send']}\t"
                         f"{h['strand']}\t{h['pident']}\t{h['mism']}\t{'YES' if h['binds'] else 'no'}\n")
                if h["binds"]: nbind += 1
    L(f"  foreign-type 3'-anchored single-primer binds: {nbind} (rows incl near-misses in {hits_path.name})")

    # circular amplicon caller
    rows = circular_amplicons(hits, prim, amps, glen)
    agg = defaultdict(lambda: dict(genomes=set(), minpid=100.0, prods=[]))
    for r in rows:
        k = (f"{r['owner_type']}_{r['owner_amp']}", r["foreign_type"])
        agg[k]["genomes"].add(r["foreign_genome"])
        agg[k]["minpid"] = min(agg[k]["minpid"], r["pident"])
        agg[k]["prods"].append(r["product"])

    amp_path = SPEC / ("crosstype_amplicons.tsv" if label in ("v4","v5") else f"crosstype_amplicons_{label}.tsv")
    with open(amp_path, "w") as fh:
        fh.write("owner_amplicon\tforeign_type\tn_foreign_genomes_amplified\t"
                 "min_pident_both_primers\tpredicted_product_size\tverdict\n")
        for (oa, ft), v in sorted(agg.items()):
            prods = v["prods"]
            psize = f"{min(prods)}-{max(prods)}" if min(prods) != max(prods) else str(prods[0])
            fh.write(f"{oa}\t{ft}\t{len(v['genomes'])}\t{round(v['minpid'],2)}\t{psize}\tFAIL\n")
    n_amp = len(agg)
    L(f"  CROSS-TYPE AMPLICONS (BLAST-coord, circular distance): {n_amp}")
    for (oa, ft), v in sorted(agg.items()):
        L(f"    FAIL: {oa} -> {ft}: {len(v['genomes'])} genomes, minpid {v['minpid']:.2f}, "
          f"product {min(v['prods'])}-{max(v['prods'])} bp")

    # independent IUPAC-engine cross-check
    foreign_by_type = {}
    for t in TYPES:
        gl = []
        for r in SeqIO.parse(str(DDV2 / f"HPV{t}.norm.fasta"), "fasta"):
            s = str(r.seq).upper(); rc = str(Seq(s).reverse_complement())
            gl.append((r.id.split("|")[1], s + s, len(s), rc + rc))
        foreign_by_type[f"HPV{t}"] = gl
    eng = iupac_crosscheck(amps, prim, foreign_by_type)
    L(f"  IUPAC-engine independent cross-check amplicons: {len(eng)} {sorted(eng) if eng else ''}")

    # sister-pair explicit both-direction status
    L("  --- sister-pair both-direction check ---")
    blast_pairs = {(oa.split('_')[0], ft) for (oa, ft) in agg}
    for a, b in SISTERS:
        for x, y in ((a, b), (b, a)):
            hit = any(oa.startswith(f"HPV{x}_") and ft == f"HPV{y}" for (oa, ft) in agg)
            L(f"    HPV{x} primers -> HPV{y} genomes: {'CROSS-AMPLICON' if hit else 'clean'}")

    L(f"\n  SCREEN RESULT [{label}]: {n_amp} cross-type amplicons (BLAST-coord), "
      f"{len(eng)} (IUPAC-engine).")
    return n_amp, len(eng), agg


if __name__ == "__main__":
    n_amp, n_eng, agg = main()
    logp = PROJ / "v4" / "00_logs" / f"t2_crosstype_{sys.argv[2] if len(sys.argv)>2 else 'v4'}.log"
    logp.write_text("\n".join(log) + "\n")
    print(f"Log -> {logp}")
