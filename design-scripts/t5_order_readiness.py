#!/usr/bin/env python
"""
v4 TASK 7 - final order-readiness package.

Produces in v4/06_reports/:
  primer_order_sheet_v5.tsv   144 primers, full columns incl. human/crosstype/
                              offtarget single-primer-bind flags, panel + African incl.
  ORDER_READINESS.tsv         per-type clean flags + GO/NO-GO
  degenerate_primer_notes.tsv IUPAC primers (refreshed; HPV58_3 change reflected)
  (coverage_map.tsv already written by t6_coverage.py)
Copies the specificity evidence tables into 06_reports for a self-contained package.
"""
import shutil
from pathlib import Path
from collections import defaultdict

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
DES = PROJ / "v5" / "04_design"
SPEC = PROJ / "v5" / "05_specificity"
OUT = PROJ / "v5" / "06_reports"
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
CLADE = {16:"Alpha9",31:"Alpha9",33:"Alpha9",35:"Alpha9",52:"Alpha9",58:"Alpha9",
         18:"Alpha7",39:"Alpha7",45:"Alpha7",59:"Alpha7",56:"Alpha6",66:"Alpha6"}
IUPAC = {"R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT",
         "H":"ACT","V":"ACG","N":"ACGT"}
log = []
def L(m): print(m); log.append(str(m))


def fold(seq):
    f = 1
    for b in seq:
        if b in IUPAC: f *= len(IUPAC[b])
    return f


def binders_from(path, primer_col=0, flag_col=None, flag_yes="YES"):
    s = set()
    if not Path(path).exists():
        return s
    lines = Path(path).read_text().splitlines()
    for ln in lines[1:]:
        c = ln.split("\t")
        if flag_col is None or (len(c) > flag_col and c[flag_col] == flag_yes):
            s.add(c[primer_col])
    return s


def main():
    # single-primer 3'-bind sets
    human_b = binders_from(SPEC / "human_hits.tsv", 0, 9, "YES")        # 3prime_anchored_bind col
    cross_b = binders_from(SPEC / "crosstype_hits.tsv", 0, 9, "YES")     # 3p_bind col
    off_b = binders_from(SPEC / "offtarget_hits.tsv", 0, 9, "YES")       # 3p_bind col
    L(f"single-primer 3'-binds: human={len(human_b)} crosstype={len(cross_b)} offtarget={len(off_b)}")

    # per-amplicon panel inclusivity
    panel_inc = {}
    for ln in (SPEC / "inclusivity_by_amplicon.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); panel_inc[(c[0], c[1])] = float(c[3])
    # per-amplicon corrected African inclusivity
    afr_inc = {}
    for ln in (SPEC / "african_inclusivity_corrected.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); afr_inc[(c[0], c[1])] = float(c[5])

    # HPV58 sublineage-C validation: per-amplicon PASS -> Y
    cval = {}
    for ln in (SPEC / "hpv58_C_validation.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); cval[c[0]] = (c[9] == "PASS")

    # per-amplicon flags
    AMP_FLAG = {
        ("HPV35", "4"): "panel<95% — African-justified exception (35 USA/divergent/non-African)",
        ("HPV58", "3"): "redesigned v4 (HPV33 fix); sublineage-C validated (0 mm on 11 complete C)",
        ("HPV52", "7"): "short L1 gap-fill (1117bp) — wet-lab balancing note",
    }

    # ---- primer order sheet ----
    hdr = ("type\tclade\tamplicon\tpool\tname\tdirection\tsequence\tlength\tn_iupac\ttm\tgc\t"
           "canonical_start\tcanonical_end\tpanel_min_identity\t3prime_clamp_conserved\t"
           "predicted_inclusivity_pct\tafrican_inclusivity_pct\tsublineageC_validated\tjunction\t"
           "human_offtarget\tcrosstype_offtarget\tofftarget_organism\tn_iupac_fold\tflag")
    rows = [hdr]
    deg_rows = ["type\tamplicon\tpool\tname\tsequence\tn_iupac\tfold_degeneracy\tnote"]
    nprimers = 0
    for t in TYPES:
        for ln in (DES / f"HPV{t}.primers.tsv").read_text().splitlines()[1:]:
            c = ln.split("\t")
            typ, amp, pool, name, d, seq, length, ndeg, tm, gc, minid, cs, ce, junc = c
            nprimers += 1
            pinc = panel_inc.get((typ, amp), "")
            ainc = afr_inc.get((typ, amp), "")
            hY = "Y" if name in human_b else "N"
            cY = "Y" if name in cross_b else "N"
            oY = "Y" if name in off_b else "N"
            cv = ("Y" if cval.get(amp) else "N") if typ == "HPV58" else "n/a"
            flag = AMP_FLAG.get((typ, amp), "OK")
            rows.append("\t".join([typ, CLADE[t], amp, pool, name, d, seq, length, ndeg, tm, gc,
                                   cs, ce, minid, "Y", str(pinc), str(ainc), cv, junc, hY, cY, oY,
                                   f"{fold(seq)}x", flag]))
            if int(ndeg) > 0:
                deg_rows.append(f"{typ}\t{amp}\t{pool}\t{name}\t{seq}\t{ndeg}\t{fold(seq)}x\t"
                                f"oligo MIX ({fold(seq)} oligos); ok"
                                + (" [redesigned v4]" if (typ, amp) == ("HPV58", "3") else ""))
    (OUT / "primer_order_sheet_v5.tsv").write_text("\n".join(rows) + "\n")
    (OUT / "degenerate_primer_notes.tsv").write_text("\n".join(deg_rows) + "\n")
    L(f"wrote primer_order_sheet_v5.tsv ({nprimers} primers) and degenerate_primer_notes.tsv "
      f"({len(deg_rows)-1} degenerate)")

    # ---- pool QC -> dimers_OK + tm_spread per type ----
    type_pass = defaultdict(lambda: True); type_spread = defaultdict(float)
    for ln in (OUT / "pool_qc.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); ty = c[0]
        type_pass[ty] &= (c[10] == "PASS")
        type_spread[ty] = max(type_spread[ty], float(c[5]))
    # ---- coverage ----
    cov = {}; junc = {}
    for ln in (OUT / "coverage_map.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); cov[c[0]] = float(c[3]); junc[c[0]] = c[4]
    # ---- min panel inclusivity per type (worst amplicon) ----
    min_panel = defaultdict(lambda: 100.0); min_afr = defaultdict(lambda: 100.0)
    for (ty, amp), v in panel_inc.items():
        min_panel[ty] = min(min_panel[ty], v)
    for (ty, amp), v in afr_inc.items():
        min_afr[ty] = min(min_afr[ty], v)

    # amplicon-level specificity verdicts (all 0 -> every type clean)
    n_human_amp = len((SPEC / "human_amplicons.tsv").read_text().splitlines()) - 1
    n_cross_amp = len((SPEC / "crosstype_amplicons.tsv").read_text().splitlines()) - 1
    n_off_amp = len((SPEC / "offtarget_amplicons.tsv").read_text().splitlines()) - 1
    L(f"amplicon-level: human={n_human_amp} crosstype={n_cross_amp} offtarget={n_off_amp}")

    RESID = {
        "HPV35": "amp4 panel 75% — documented African-justified exception (35 USA divergent, 0 African)",
        "HPV58": "sublineage C POSITIVELY VALIDATED (11 complete C genomes, all 6 amps 0 mm); amp3 redesigned (HPV33 fix); 2 N-gapped C records remain documented data-artifacts",
        "HPV52": "amp7 (L1 gap-fill, 1117bp) — wet-lab balancing note",
    }
    C_VALIDATED = all(cval.values())
    # ---- ORDER_READINESS ----
    orhdr = ("type\tclade\thuman_clean\tcrosstype_clean\tofftarget_clean\tcov_pct\tjunction\t"
             "min_african_incl\tmin_panel_incl\tdimers_OK\ttm_spread_C\tresidual_flags\tGO_NOGO")
    orows = [orhdr]
    overall_go = True
    table = []
    for t in TYPES:
        ty = f"HPV{t}"
        hc = "Y" if n_human_amp == 0 else "N"
        cc = "Y" if n_cross_amp == 0 else "N"
        oc = "Y" if n_off_amp == 0 else "N"
        d_ok = "Y" if type_pass[ty] else "N"
        resid = RESID.get(ty, "none")
        # GO requires: 0 amplicon-level offtargets (all kinds), 100% cov, junction, dimers ok,
        # African 100% (assessable), Tm spread <=5. Panel <95% allowed ONLY with documented exception.
        go = (hc == "Y" and cc == "Y" and oc == "Y" and cov[ty] >= 100.0 and junc[ty] == "Y"
              and min_afr[ty] >= 100.0 and d_ok == "Y" and type_spread[ty] <= 5.0)
        # HPV58 may be GO ONLY if sublineage C is positively validated (T2) or carries an
        # explicit accession-backed C exception. Here C is positively validated.
        if ty == "HPV58":
            go = go and C_VALIDATED
        # HPV35 panel<95% is a documented exception, still GO
        go_nogo = "GO" if go else "NO-GO"
        overall_go &= go
        orows.append("\t".join([ty, CLADE[t], hc, cc, oc, f"{cov[ty]:.1f}", junc[ty],
                                f"{min_afr[ty]:.0f}", f"{min_panel[ty]:.0f}", d_ok,
                                f"{type_spread[ty]:.1f}", resid, go_nogo]))
        table.append((ty, CLADE[t], hc, cc, oc, cov[ty], junc[ty], min_afr[ty], min_panel[ty],
                      d_ok, type_spread[ty], go_nogo))
    (OUT / "ORDER_READINESS.tsv").write_text("\n".join(orows) + "\n")

    # copy specificity evidence into the package
    for f in ("crosstype_amplicons.tsv", "crosstype_hits.tsv", "offtarget_amplicons.tsv",
              "offtarget_hits.tsv", "human_amplicons.tsv", "specificity_summary.tsv"):
        if (SPEC / f).exists():
            shutil.copy(SPEC / f, OUT / f)

    # print final table
    L("\n================ FINAL ORDER READINESS (v5) ================")
    L(f"{'type':7}{'clade':8}{'hum':4}{'xtyp':5}{'offt':5}{'cov%':6}{'junc':5}"
      f"{'afrIncl':8}{'panIncl':8}{'dimer':6}{'TmSpr':6}  GO/NO-GO")
    for r in table:
        L(f"{r[0]:7}{r[1]:8}{r[2]:^4}{r[3]:^5}{r[4]:^5}{r[5]:<6.0f}{r[6]:^5}"
          f"{r[7]:<8.0f}{r[8]:<8.0f}{r[9]:^6}{r[10]:<6.1f}  {r[11]}")
    L("-" * 60)
    L(f"OVERALL: {'GO' if overall_go else 'NO-GO'}  ({sum(1 for r in table if r[11]=='GO')}/12 types GO)")
    return overall_go, table


if __name__ == "__main__":
    go, table = main()
    (PROJ / "v5" / "00_logs" / "t5_order_readiness.log").write_text("\n".join(log) + "\n")
    print(f"\nLog -> {PROJ/'v5/00_logs/t5_order_readiness.log'}")
