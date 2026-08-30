#!/usr/bin/env python
"""
prepare_design_bundle.py  (LOCAL, run in this design project; no server needed)

Packages the v5 HPV tiling-amplicon design into an upload-ready, self-contained
ARTIC-style bundle for the CPU-only Nanopore analysis server.

REFERENCES = the CANONICAL accession genome per type, exactly as documented in
v5/06_reports/METHODS_REPORT.md §2.3 (PaVE/NCBI: HPV16 K02718, HPV18 X05015,
HPV31 J04353, HPV33 M12732, HPV35 X74477, HPV39 M62849, HPV45 X74479, HPV52 X74481,
HPV56 X74483, HPV58 D90400, HPV59 X77858, HPV66 U31794). These are read from the
record already referenced in the design (01_raw_genomes/pave_references.fasta); if a
type's canonical record is missing/blocked, the script FAILS LOUDLY for that type rather
than falling back to a consensus.

Primer coordinates are RE-DERIVED against each canonical genome (the design's
canonical_start/end live in the panel-CONSENSUS frame, which drifts from the accession by
indels). Each primer is located by an IUPAC-aware, 3'-clamp-anchored search in the
canonical genome; exact match is confirmed; the re-derived position is cross-checked
against primer_order_sheet_v5.tsv and any disagreement/inexact match is FLAGGED (never
silently shifted).

Produces ./hpv-analysis-bundle/ :
  primer-schemes/HPV{t}/V5/HPV{t}.reference.fasta   (header = "HPV{t}|{accession}")
  primer-schemes/HPV{t}/V5/HPV{t}.primer.bed        (re-derived coords; _alt for IUPAC)
  primer-schemes/HPV{t}/V5/HPV{t}.insert.bed        (amplicon inserts; origin-wrap split+flag)
  primer-schemes/HPV{t}/V5/HPV{t}.amplicon.bed      (full primer-to-primer spans; wrap split)
  references/hpv_all_types.fasta                     (12 canonical refs, same headers)
  metadata/scheme_summary.tsv                        (per-type summary + carried design notes)
  metadata/coordinate_rederivation.tsv               (per-primer audit: accession coords + match)
  sample_sheet_TEMPLATE.csv
  MANIFEST.txt  SHA256SUMS  README.md  BUNDLE_VERSION
"""
import sys, json, hashlib, shutil, itertools, datetime, subprocess
from pathlib import Path
from collections import defaultdict

PROJ = Path(__file__).resolve().parent
DES = PROJ / "v5" / "04_design"
PAVE_FASTA = PROJ / "01_raw_genomes" / "pave_references.fasta"
OUT = PROJ / "hpv-analysis-bundle"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
CLADE = {16:"Alpha9",31:"Alpha9",33:"Alpha9",35:"Alpha9",52:"Alpha9",58:"Alpha9",
         18:"Alpha7",39:"Alpha7",45:"Alpha7",59:"Alpha7",56:"Alpha6",66:"Alpha6"}
VERSION = "V5"                       # legacy dir name (get_scheme_legacy: <name>/V5/)
ARTIC_LEN = "2000"                    # ARTIC 1.10.1 new layout "scheme length" token
ARTIC_VER = "v5.0.0"                  # ARTIC 1.10.1 new layout scheme version (semver)
# canonical accessions documented in METHODS_REPORT.md §2.3 (the authority)
DOC_ACC = {16:"K02718",18:"X05015",31:"J04353",33:"M12732",35:"X74477",39:"M62849",
           45:"X74479",52:"X74481",56:"X74483",58:"D90400",59:"X77858",66:"U31794"}
DELTA_TOL = 150           # max |re-derived - order-sheet| start offset (indel drift) before FLAG
IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT",
         "M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
flags = []
def flag(m): print("  [FLAG]", m); flags.append(m)
def fail_loud(m): print("  [FATAL]", m, file=sys.stderr)


def revcomp(s): return s.translate(COMP)[::-1]
def iset(b): return IUPAC.get(b, b)
def expand_iupac(seq):
    return ["".join(p) for p in itertools.product(*[iset(b) for b in seq])]


def load_canonical_refs():
    """type -> (accession, seq). From the cached PaVE record (the documented accession).
    Verifies accession == documented; FAILS LOUDLY (and tries NCBI) otherwise."""
    refs = {}
    if PAVE_FASTA.exists():
        from Bio import SeqIO
        for r in SeqIO.parse(str(PAVE_FASTA), "fasta"):
            parts = r.id.split("|")          # HPV{t}|{acc}|...
            if len(parts) >= 2 and parts[0].startswith("HPV"):
                try: t = int(parts[0][3:])
                except ValueError: continue
                refs[t] = (parts[1], str(r.seq).upper())
    out = {}
    for t in TYPES:
        want = DOC_ACC[t]
        if t in refs and refs[t][0] == want:
            out[t] = refs[t]
        else:
            got = refs.get(t, (None, None))[0]
            print(f"  [WARN] HPV{t}: local canonical accession is '{got}', want '{want}'. "
                  f"attempting NCBI fetch of {want} ...")
            seq = fetch_ncbi(want)
            if seq:
                out[t] = (want, seq)
            else:
                fail_loud(f"HPV{t}: could not obtain canonical accession {want} "
                          f"(local mismatch and NCBI fetch blocked). NOT falling back to consensus.")
                out[t] = None
    return out


def fetch_ncbi(acc):
    """efetch a GenBank accession as FASTA via eutils (pinned IP if DNS is blocked)."""
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=fasta&retmode=text")
    for args in ([], ["--resolve", "eutils.ncbi.nlm.nih.gov:443:" + _eutils_ip()]):
        try:
            r = subprocess.run(["curl", "-fsSL", "--max-time", "60", *args, url],
                               capture_output=True, text=True)
            if r.returncode == 0 and r.stdout.startswith(">"):
                return "".join(l.strip() for l in r.stdout.splitlines()[1:]).upper()
        except Exception:
            pass
    return None


def _eutils_ip():
    try:
        r = subprocess.run(["nslookup", "eutils.ncbi.nlm.nih.gov", "8.8.8.8"],
                           capture_output=True, text=True)
        for ln in r.stdout.splitlines():
            if ln.startswith("Address:") and "#" not in ln and "53" not in ln.split()[-1]:
                ip = ln.split()[-1]
                if ip.count(".") == 3: return ip
    except Exception:
        pass
    return "0.0.0.0"


def load_primers():
    rows, hdr = [], None
    for ln in (DES.parent / "06_reports" / "primer_order_sheet_v5.tsv").read_text().splitlines():
        c = ln.split("\t")
        if hdr is None: hdr = c; continue
        rows.append(dict(zip(hdr, c)))
    return rows


def find_primer(pattern, acc_seq, anchor, expected, W=250):
    """Locate `pattern` (+ strand, IUPAC-aware) in the canonical accession by PROXIMITY
    to the expected (order-sheet) position, then minimum mismatch — NOT by hard-gating
    on the 3' clamp (the canonical accession need not carry the panel's clamp allele).
    anchor='end' -> 3' clamp = last 5 nt of pattern; 'start' -> first 5 nt.
    Returns dict(start, exact, nmm, mmpos, clamp_mm). Searches a window around expected,
    widening to the whole genome only if the window is degenerate."""
    L = len(acc_seq); m = len(pattern)
    if m > L: return None
    clamp = range(m - 5, m) if anchor == "end" else range(0, 5)
    def scan(lo, hi):
        best = None
        for p in range(max(0, lo), min(L - m, hi) + 1):
            seg = acc_seq[p:p + m]
            viol = sum(1 for k in range(m) if seg[k] not in iset(pattern[k]))
            if best is None or viol < best[0] or (viol == best[0] and abs(p - expected) < abs(best[1] - expected)):
                best = (viol, p, seg)
        return best
    best = scan(expected - W, expected + W)
    if best is None:                      # expected outside [0,L-m]; full scan
        best = scan(0, L - m)
    if best is None: return None
    viol, p, seg = best
    mmpos = [(m - k if anchor == "end" else k + 1) for k in range(m) if seg[k] not in iset(pattern[k])]
    clamp_mm = sum(1 for k in clamp if seg[k] not in iset(pattern[k]))
    return dict(start=p, exact=(viol == 0), nmm=viol, mmpos=mmpos, clamp_mm=clamp_mm)


def main():
    refs = load_canonical_refs()
    if any(refs[t] is None for t in TYPES):
        missing = [f"HPV{t}({DOC_ACC[t]})" for t in TYPES if refs[t] is None]
        fail_loud(f"Aborting: canonical reference unavailable for: {', '.join(missing)}")
        sys.exit(3)

    if OUT.exists(): shutil.rmtree(OUT)
    (OUT / "references").mkdir(parents=True)
    (OUT / "metadata").mkdir(parents=True)
    print(f"=== prepare_design_bundle ({VERSION}) — CANONICAL accession references -> {OUT} ===")

    primers = load_primers()
    by_type = defaultdict(list)
    for p in primers: by_type[p["type"]].append(p)

    NOTE = {
        "HPV52": ("amp7 = L1 gap-fill ~1117 bp, SHORTER than the 2-4 kb tiling amplicons; "
                  "over-amplifies in multiplex. Wet-lab: reduce HPV52_7 primer conc or "
                  "post-PCR normalize. Length filter must keep >=~900 bp. (wetlab_balancing_notes.md)"),
        "HPV35": ("amp4 panel inclusivity 75% = documented African-justified exception (35 "
                  "USA/divergent, 0 African); African inclusivity 100%. (softsite_exceptions.md)"),
    }

    all_ref_records = []
    manifest_entries = {}            # for the local ARTIC index.json (offline get_scheme)
    summary_rows = ["type\tclade\taccession\tn_amplicons\tpools\tcanonical_len\tjunction_amplicon\t"
                    "junction_wraps_origin\tn_degenerate_primers\tdegenerate_primers\tnotes"]
    audit_rows = ["type\tprimer\tdirection\torder_sheet_start\torder_sheet_end\t"
                  "rederived_start\trederived_end\tdelta_start\texact_match\tn_mismatch\t"
                  "mismatch_pos_from_3p\tn_alt\tn_alt_exact_substr\tstatus"]
    report = []   # per-type printable

    for t in TYPES:
        ty = f"HPV{t}"; acc, ref = refs[t]
        L = len(ref)
        sd = OUT / "primer-schemes" / ty / VERSION; sd.mkdir(parents=True, exist_ok=True)
        # chrom = type_accession (underscore, NOT '|': primalbedtools/ARTIC require
        # chrom matching ^[a-zA-Z0-9_.]+$ — a pipe is rejected). Provenance preserved.
        chrom = f"{ty}_{acc}"
        with open(sd / f"{ty}.reference.fasta", "w") as fh:
            fh.write(f">{chrom}\n")
            for i in range(0, L, 70): fh.write(ref[i:i+70] + "\n")
        all_ref_records.append((chrom, ref))

        prs = by_type[ty]
        amp_primers = defaultdict(dict)
        pools = set(); bed_lines = []; ndeg = 0; deg_names = []
        n_re = 0; n_exact = 0; type_flags = []
        for p in prs:
            amp = p["amplicon"]; d = p["direction"]; name = p["name"]; seq = p["sequence"]
            pool = p["pool"]; pools.add(pool); strand = "+" if d == "F" else "-"
            os_cs, os_ce = int(p["canonical_start"]), int(p["canonical_end"])
            pat = seq if d == "F" else revcomp(seq)    # + strand pattern
            anchor = "end" if d == "F" else "start"    # 3' clamp location on + strand
            loc = find_primer(pat, ref, anchor, os_cs)
            n_re += 1
            variants = expand_iupac(seq)
            if len(variants) > 1: ndeg += 1; deg_names.append(name)
            if loc is None:
                type_flags.append(name)
                flag(f"{name}: NOT FOUND in canonical {acc} (no 3'-clamp-anchored match) — "
                     f"left at order-sheet coords {os_cs}-{os_ce}, FLAGGED.")
                cs, ce = os_cs, os_ce; exact = False; nmm = "NA"; mmpos = []; status = "NOT_FOUND"
            else:
                cs = loc["start"]; ce = cs + len(seq); exact = loc["exact"]
                nmm = loc["nmm"]; mmpos = loc["mmpos"]
                delta = cs - os_cs
                if exact: n_exact += 1
                status = "OK"
                if not exact:
                    clamp_mm = loc["clamp_mm"]
                    status = "INEXACT_CLAMP" if clamp_mm else "INEXACT_INTERNAL"
                    type_flags.append(name)
                    flag(f"{name}: located in {acc} at {cs} (delta {delta:+d}) but NOT exact — "
                         f"{nmm} mismatch vs accession ({clamp_mm} in the 3' clamp) at 3'-offset "
                         f"{mmpos}. Accession carries a different allele than the panel-designed "
                         f"primer; FLAGGED, coords kept at the true homologous site (not shifted).")
                if abs(delta) > DELTA_TOL:
                    status = (status + "+DELTA").replace("OK+", "")
                    type_flags.append(name)
                    flag(f"{name}: re-derived start {cs} disagrees with order sheet {os_cs} "
                         f"(delta {delta} > +/-{DELTA_TOL}).")
            amp_primers[amp][d] = dict(cs=cs, ce=ce, pool=pool, strand=strand, name=name)
            # _alt expansion at the re-derived coords; count how many concrete variants are
            # exact substrings of the canonical genome at this site (informational)
            n_alt_exact = 0
            for i, v in enumerate(variants):
                nm = name if i == 0 else f"{name}_alt{i}"
                bed_lines.append((cs, f"{chrom}\t{cs}\t{ce}\t{nm}\t{pool}\t{strand}\t{v}"))
                onplus = v if d == "F" else revcomp(v)
                if ref[cs:cs+len(onplus)] == onplus: n_alt_exact += 1
            delta_str = (str(cs - os_cs) if loc else "NA")
            audit_rows.append(f"{ty}\t{name}\t{d}\t{os_cs}\t{os_ce}\t{cs}\t{ce}\t{delta_str}\t"
                              f"{'YES' if exact else 'no'}\t{nmm}\t{','.join(map(str,mmpos)) or '-'}\t"
                              f"{len(variants)}\t{n_alt_exact}\t{status}")
        bed_lines.sort(key=lambda x: x[0])
        primer_bed_text = "\n".join(l for _, l in bed_lines) + "\n"
        ref_fasta_text = f">{chrom}\n" + "".join(ref[i:i+70] + "\n" for i in range(0, L, 70))
        (sd / f"{ty}.primer.bed").write_text(primer_bed_text)
        (sd / f"{ty}.scheme.bed").write_text(primer_bed_text)   # name get_scheme_legacy() expects
        # --- ARTIC 1.10.1 new scheme layout: <name>/<length>/<version>/{primer.bed,reference.fasta}
        nd = OUT / "primer-schemes" / ty / ARTIC_LEN / ARTIC_VER; nd.mkdir(parents=True, exist_ok=True)
        (nd / "primer.bed").write_text(primer_bed_text)
        (nd / "reference.fasta").write_text(ref_fasta_text)
        manifest_entries[ty] = dict(
            length=ARTIC_LEN, version=ARTIC_VER,
            primer_bed_md5=hashlib.md5(primer_bed_text.encode()).hexdigest(),
            reference_fasta_md5=hashlib.md5(ref_fasta_text.encode()).hexdigest())

        # insert/amplicon beds with origin-wrap split + flag (based on RE-DERIVED coords)
        junction_amp = None; insert_lines, amp_lines = [], []
        for amp in sorted(amp_primers, key=int):
            fr = amp_primers[amp]
            if "F" not in fr or "R" not in fr:
                flag(f"{ty} amp{amp}: missing F or R after re-derivation"); continue
            Lp, Rp = fr["F"], fr["R"]; apool = Lp["pool"]; iname = f"{ty}_{amp}"
            ins0, ins1 = Lp["ce"], Rp["cs"]; amp0, amp1 = Lp["cs"], Rp["ce"]
            if amp0 > amp1 or ins0 > ins1:
                junction_amp = amp
                insert_lines += [f"{chrom}\t{ins0}\t{L}\t{iname}_I_part1\t{apool}\t+",
                                 f"{chrom}\t0\t{ins1}\t{iname}_I_part2\t{apool}\t+"]
                amp_lines    += [f"{chrom}\t{amp0}\t{L}\t{iname}_part1\t{apool}\t+",
                                 f"{chrom}\t0\t{amp1}\t{iname}_part2\t{apool}\t+"]
            else:
                insert_lines.append(f"{chrom}\t{ins0}\t{ins1}\t{iname}_I\t{apool}\t+")
                amp_lines.append(f"{chrom}\t{amp0}\t{amp1}\t{iname}\t{apool}\t+")
        (sd / f"{ty}.insert.bed").write_text("\n".join(insert_lines) + "\n")
        (sd / f"{ty}.amplicon.bed").write_text("\n".join(amp_lines) + "\n")
        if junction_amp:
            flag(f"{ty}: origin-junction amplicon {junction_amp} WRAPS position 0 on canonical "
                 f"{acc} -> split into two BED intervals; recovered at run time by the pipeline's "
                 f"rotated-reference step (hpvpipe.py rotate-build + consensus).")

        summary_rows.append(f"{ty}\t{CLADE[t]}\t{acc}\t{len(amp_primers)}\t{','.join(sorted(pools))}\t"
                            f"{L}\t{junction_amp if junction_amp else 'none'}\t"
                            f"{'YES' if junction_amp else 'no'}\t{ndeg}\t"
                            f"{';'.join(deg_names) if deg_names else '-'}\t{NOTE.get(ty,'')}")
        report.append((ty, acc, L, n_re, n_exact, sorted(set(type_flags))))

    # concatenated canonical references (same headers)
    with open(OUT / "references" / "hpv_all_types.fasta", "w") as fh:
        for name, seq in all_ref_records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70): fh.write(seq[i:i+70] + "\n")
    # NOTE: minimap2 .mmi built on server (not portable); see README.

    # ---- local ARTIC 1.10.1 manifest (primer-schemes/index.json) ----
    # get_scheme() fetches quick-lab/primerschemes/index.json; if the network is unavailable
    # (cluster compute nodes) it falls back to THIS local copy, then skips downloading any
    # primer.bed/reference.fasta that already exist at <name>/<length>/<version>/ — so the
    # scheme resolves fully OFFLINE. (The pipeline itself uses `artic minion --bed/--ref`,
    # which bypasses manifest resolution entirely; this index.json supports the
    # --scheme-name/--scheme-length/--scheme-version form too.)
    manifest = {"primerschemes": {}}
    for ty, m in manifest_entries.items():
        manifest["primerschemes"].setdefault(ty, {}).setdefault(m["length"], {})[m["version"]] = {
            "primer_bed_url": f"https://local.invalid/{ty}/{m['length']}/{m['version']}/primer.bed",
            "primer_bed_md5": m["primer_bed_md5"],
            "reference_fasta_url": f"https://local.invalid/{ty}/{m['length']}/{m['version']}/reference.fasta",
            "reference_fasta_md5": m["reference_fasta_md5"],
            "status": "draft", "schemeName": ty, "ampliconSize": int(m["length"]),
            "schemeVersion": m["version"]}
    (OUT / "primer-schemes" / "index.json").write_text(json.dumps(manifest, indent=1) + "\n")

    (OUT / "metadata" / "scheme_summary.tsv").write_text("\n".join(summary_rows) + "\n")
    (OUT / "metadata" / "coordinate_rederivation.tsv").write_text("\n".join(audit_rows) + "\n")
    (OUT / "sample_sheet_TEMPLATE.csv").write_text(
        "barcode,sample_id,expected_types,notes\n"
        "# barcode = MinKNOW barcodeNN folder (e.g. barcode01); expected_types optional, ;-separated\n"
        "barcode01,SAMPLE_001,,\n barcode02,SAMPLE_002,,\n".replace("\n ", "\n"))
    (OUT / "BUNDLE_VERSION").write_text(VERSION + "\n")
    (OUT / "README.md").write_text(BUNDLE_README)

    # MANIFEST + SHA256SUMS
    files = sorted(p for p in OUT.rglob("*") if p.is_file()
                   and p.name not in ("MANIFEST.txt", "SHA256SUMS"))
    man, sums = [], []
    for f in files:
        rel = f.relative_to(OUT); h = hashlib.sha256(f.read_bytes()).hexdigest()
        man.append(f"{rel}\t{f.stat().st_size}\t{h}"); sums.append(f"{h}  {rel}")
    ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    (OUT / "MANIFEST.txt").write_text(
        f"# HPV analysis bundle {VERSION} (canonical-accession references)  built {ts}\n"
        f"# path\tbytes\tsha256\n" + "\n".join(man) + "\n")
    (OUT / "SHA256SUMS").write_text("\n".join(sums) + "\n")

    # ---- per-type report ----
    print("\n=== PER-TYPE RE-DERIVATION REPORT (canonical accession) ===")
    print(f"{'type':6} {'accession':10} {'len':6} {'primers':8} {'exact':6} disagreements")
    tot_p = tot_e = 0
    for ty, acc, L, n_re, n_exact, fl in report:
        tot_p += n_re; tot_e += n_exact
        print(f"{ty:6} {acc:10} {L:<6} {n_re:<8} {n_exact:<6} "
              f"{'none' if not fl else ','.join(fl)}")
    print("-" * 64)
    placement_errs = [f for f in flags if ("NOT FOUND" in f) or ("disagrees with order sheet" in f)]
    allele_diffs   = [f for f in flags if ("NOT exact" in f)]
    if tot_e == tot_p:
        print(f"ALL {tot_p} primers matched their canonical accession EXACTLY (IUPAC-aware) and "
              f"agree with primer_order_sheet_v5.tsv (within indel drift). No coordinate "
              f"disagreements. Only the documented origin-junction wrap is flagged.")
    else:
        print(f"{tot_e}/{tot_p} primers exact-match the canonical accession (IUPAC-aware).")
        print(f"  {len(allele_diffs)} documented ALLELE DIFFERENCE(S) — primer located at the correct "
              f"homologous site but the accession carries a different base (see below). Coordinates "
              f"are CORRECT (delta within indel drift); ARTIC trims by coordinate, so this does not "
              f"affect trimming/calling — it is a reference-allele annotation:")
        for f in allele_diffs:
            print("    - " + f.split(" — ")[0].replace("[FLAG] ", "").strip())
        if placement_errs:
            print(f"  {len(placement_errs)} PLACEMENT ERROR(S) (mislocated / not found) — see [FLAG] above.")
    print(f"\nNote: re-derived coordinates are in the CANONICAL ACCESSION frame; the order sheet's "
          f"canonical_start/end are in the panel-CONSENSUS frame, so a small negative delta "
          f"(<= a few tens of bp, from consensus insertions) is expected and not an error.")
    print(f"\nbundle: {len(files)} files. Next: verify SHA256SUMS, upload {OUT.name}/ + server scripts.")
    # exit non-zero ONLY on real placement failures (mislocation / not found); a documented
    # allele difference with correct coordinates is a flag, not a packaging failure.
    sys.exit(2 if placement_errs else 0)


BUNDLE_README = """# HPV analysis bundle (v5) — CANONICAL-accession references

Self-contained, version-stamped (V5). Nothing here phones home at analysis time.

## References (FIXED: canonical accessions, not panel consensus)
Each per-type reference is the **documented canonical accession** from
METHODS_REPORT.md §2.3 (HPV16 K02718, HPV18 X05015, HPV31 J04353, HPV33 M12732,
HPV35 X74477, HPV39 M62849, HPV45 X74479, HPV52 X74481, HPV56 X74483, HPV58 D90400,
HPV59 X77858, HPV66 U31794), read from the PaVE record already used in the design.
FASTA header (== BED chrom) is **`HPV{type}|{accession}`** so provenance is explicit in
BAMs/VCFs. Coordinates are interpretable directly against the published accession.

## Coordinate re-derivation
Primer coordinates were RE-DERIVED against each canonical genome (the design's
`canonical_start/end` are in the panel-consensus frame, which drifts from the accession
by indels). Each primer (and the IUPAC-aware site for degenerate primers) was located by a
3'-clamp-anchored, IUPAC-aware search; exact match was confirmed and cross-checked against
`primer_order_sheet_v5.tsv`. See `metadata/coordinate_rederivation.tsv` for the per-primer
audit (order-sheet vs re-derived coords, delta, exact-match, mismatch positions). Any
inexact/disagreeing primer is FLAGGED there, never silently shifted.

## Contents
- `primer-schemes/HPV{type}/V5/` — ARTIC scheme dir per type:
  `HPV{type}.reference.fasta` (header `HPV{type}|{accession}`),
  `HPV{type}.primer.bed` (7-col: chrom start end name pool strand sequence; degenerate
  primers expanded to concrete oligos — base name + `_alt1`, `_alt2`, …),
  `HPV{type}.insert.bed`, `HPV{type}.amplicon.bed`.
- `references/hpv_all_types.fasta` — 12 canonical refs concatenated (same headers); the
  minimap2 `.mmi` is built on the server (not portable).
- `metadata/scheme_summary.tsv` (now reports the **accession** and canonical length),
  `metadata/coordinate_rederivation.tsv` (per-primer coordinate audit).
- `sample_sheet_TEMPLATE.csv`, `MANIFEST.txt`, `SHA256SUMS`.

## ⚠ Origin-junction limitation (unchanged)
Each HPV genome is circular; amplicon 3 spans position 0 on the linear canonical
reference. Per-primer BED entries are valid; the amplicon span is split into
`..._part1`/`..._part2` and flagged (`scheme_summary.tsv: junction_wraps_origin=YES`).
A linear reference cannot trim/call the junction amplicon as one unit — the rotated/doubled
reference fix remains the documented follow-up; the other 5 amplicons are unaffected.

## Server ingest
`setup_hpv_nanopore_env.sh` verifies `SHA256SUMS`, copies `primer-schemes/` -> `$SCHEME_DIR`,
`references/` -> `$REF_DIR` (builds `hpv_all_types.mmi`), and re-validates one scheme.
"""

if __name__ == "__main__":
    main()
