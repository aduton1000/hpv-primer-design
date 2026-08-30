#!/usr/bin/env python
"""
v5 METHODS_REPORT revision - TASK 1 (derivation only; no re-fetch/re-analysis).

Per-type x origin breakdown of the final v5 scheme genomes, using the CORRECTED
geography (country recovered from ISO-3166 alpha-3 isolate-code prefixes where the
GenBank country field is empty) - the same logic that produced the corrected African
denominators in v4/05_specificity/african_inclusivity_corrected.tsv (report S7.4).

Sources (read-only):
  v2/02_data/raw/HPV{t}.provenance.tsv   (via v4/scripts/t4b_geography_audit.geography)
  v2/02_data/sequence_provenance.tsv     (full panel n per type, status)
  v5/02_data/hpv58_C_provenance.tsv      (11 complete C genomes)
  v4/05_specificity/african_inclusivity_corrected.tsv (cross-check African n)

Outputs (tables for the report; printed + saved):
  v5/06_reports/origin_breakdown_region.tsv
  v5/06_reports/origin_breakdown_african_by_country.tsv
"""
import importlib.util, re
from pathlib import Path
from collections import defaultdict, Counter

PROJ = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
s = importlib.util.spec_from_file_location("t4b", str(PROJ / "v4/scripts/t4b_geography_audit.py"))
t4b = importlib.util.module_from_spec(s); s.loader.exec_module(t4b)
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
REGIONS = ["Africa", "Europe", "Asia", "Americas", "Oceania", "Unknown"]
log = []
def L(m): print(m); log.append(str(m))

# country (normalized) -> region. Only countries that actually appear in the files.
COUNTRY_REGION = {
    # Africa
    "Togo": "Africa", "South Africa": "Africa", "Nigeria": "Africa", "Rwanda": "Africa",
    "Guinea": "Africa", "Zimbabwe": "Africa", "Mali": "Africa", "Cameroon": "Africa",
    "Ethiopia": "Africa", "Kenya": "Africa", "Mozambique": "Africa", "Senegal": "Africa",
    "Tanzania": "Africa", "Uganda": "Africa", "Zambia": "Africa",
    # Europe
    "Netherlands": "Europe", "Luxembourg": "Europe", "Latvia": "Europe",
    # Asia
    "Japan": "Asia", "China": "Asia", "India": "Asia", "Bangladesh": "Asia",
    "South Korea": "Asia", "Thailand": "Asia",
    # Americas
    "USA": "Americas", "Mexico": "Americas", "Brazil": "Americas",
    "Argentina": "Americas", "Honduras": "Americas",
}

def norm_country(c):
    """Normalize a raw country string: drop sub-national detail after ':' , trim."""
    if not c:
        return ""
    base = c.split(":")[0].strip()
    return base

def classify(country, is_afr):
    """Return (region, normalized_country). is_afr from corrected geography."""
    nc = norm_country(country)
    if is_afr:
        return "Africa", (nc if nc else "Africa (unspecified)")
    if nc and nc in COUNTRY_REGION:
        return COUNTRY_REGION[nc], nc
    if nc:
        return "Unknown", nc          # a country we did not map -> flagged Unknown, not guessed
    return "Unknown", ""

def main():
    # full-panel n per type (sequence_provenance) for reconciliation
    full_n = Counter()
    status_by_acc = {}
    for ln in (PROJ / "v2/02_data/sequence_provenance.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); full_n[c[1]] += 1; status_by_acc[c[0]] = c[9]
    # corrected African n per type (cross-check)
    afr_corr = {}
    for ln in (PROJ / "v4/05_specificity/african_inclusivity_corrected.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t"); afr_corr.setdefault(c[0], int(c[2]))

    region_tbl = {}      # type -> Counter(region)
    afr_country = defaultdict(Counter)   # type -> Counter(african country)
    unmapped_examples = Counter()
    src_counter = Counter()
    afr_check = {}
    for t in TYPES:
        geo = t4b.geography(t)   # acc -> (country_str, is_african, source)
        rc = Counter()
        for acc, (country, is_afr, source) in geo.items():
            src_counter[source] += 1
            region, nc = classify(country, is_afr)
            rc[region] += 1
            if region == "Africa":
                afr_country[f"HPV{t}"][nc] += 1
            if region == "Unknown" and nc:
                unmapped_examples[nc] += 1
        region_tbl[f"HPV{t}"] = rc
        afr_check[f"HPV{t}"] = rc["Africa"]

    # ---- reconciliation ----
    L("=== RECONCILIATION ===")
    all_ok = True
    for t in TYPES:
        ty = f"HPV{t}"
        tot = sum(region_tbl[ty].values())
        ok = (tot == full_n[ty])
        all_ok &= ok
        ac = afr_check[ty]; ace = afr_corr.get(ty)
        afr_ok = (ac == ace)
        all_ok &= afr_ok
        L(f"  {ty}: region-table total={tot} vs full-panel n={full_n[ty]} {'OK' if ok else 'MISMATCH'}"
          f" | African={ac} vs corrected-file={ace} {'OK' if afr_ok else 'MISMATCH'}")
    grand = sum(sum(region_tbl[f'HPV{t}'].values()) for t in TYPES)
    L(f"  v5 genome TOTAL (full panel) = {grand}  (expected 1044)")
    # C-genome reconciliation
    cgen = [ln.split("\t")[0] for ln in (PROJ/"v5/02_data/hpv58_C_provenance.tsv").read_text().splitlines()[1:]
            if ln.split("\t")[-1] == "YES"]
    in_panel = sum(1 for a in cgen if a in status_by_acc)
    L(f"  HPV58 complete-C genomes: {len(cgen)}; already in the 1044 panel: {in_panel}; "
      f"new additions: {len(cgen)-in_panel}")
    L(f"  corrected-geography source mix: {dict(src_counter)}")
    L(f"  Unknown-binned country strings not in region map: {dict(unmapped_examples)}")

    # ---- TABLE (a): per-type x region ----
    L("\n=== TABLE (a) Per-type x geographic region (corrected geography) ===")
    hdr = ["type"] + REGIONS + ["Type_total", "pct_African"]
    rows_a = ["\t".join(hdr)]
    L("\t".join(f"{h:>9}" if h != "type" else f"{h:<6}" for h in hdr))
    tot_region = Counter()
    for t in TYPES:
        ty = f"HPV{t}"; rc = region_tbl[ty]
        tot = sum(rc.values())
        cells = [rc.get(r, 0) for r in REGIONS]
        pct = round(100.0 * rc.get("Africa", 0) / tot, 1) if tot else 0.0
        rows_a.append("\t".join([ty] + [str(x) for x in cells] + [str(tot), str(pct)]))
        L(f"{ty:<6}" + "".join(f"{x:>9}" for x in cells) + f"{tot:>9}" + f"{pct:>9}")
        for r in REGIONS: tot_region[r] += rc.get(r, 0)
    gtot = sum(tot_region.values())
    gpct = round(100.0 * tot_region["Africa"] / gtot, 1)
    rows_a.append("\t".join(["TOTAL"] + [str(tot_region[r]) for r in REGIONS] + [str(gtot), str(gpct)]))
    L(f"{'TOTAL':<6}" + "".join(f"{tot_region[r]:>9}" for r in REGIONS) + f"{gtot:>9}" + f"{gpct:>9}")
    (PROJ / "v5/06_reports/origin_breakdown_region.tsv").write_text("\n".join(rows_a) + "\n")

    # ---- TABLE (b): per-type African by country ----
    L("\n=== TABLE (b) African genomes by country (corrected geography) ===")
    afr_countries = sorted({c for ty in afr_country for c in afr_country[ty]})
    rows_b = ["\t".join(["type"] + afr_countries + ["African_total"])]
    L("type\t" + "\t".join(afr_countries) + "\tAfrican_total")
    col_tot = Counter()
    for t in TYPES:
        ty = f"HPV{t}"; cc = afr_country[ty]
        if sum(cc.values()) == 0:
            continue
        cells = [cc.get(c, 0) for c in afr_countries]
        rows_b.append("\t".join([ty] + [str(x) for x in cells] + [str(sum(cc.values()))]))
        L(f"{ty}\t" + "\t".join(str(x) for x in cells) + f"\t{sum(cc.values())}")
        for c in afr_countries: col_tot[c] += cc.get(c, 0)
    rows_b.append("\t".join(["TOTAL"] + [str(col_tot[c]) for c in afr_countries] + [str(sum(col_tot.values()))]))
    L("TOTAL\t" + "\t".join(str(col_tot[c]) for c in afr_countries) + f"\t{sum(col_tot.values())}")
    (PROJ / "v5/06_reports/origin_breakdown_african_by_country.tsv").write_text("\n".join(rows_b) + "\n")

    L(f"\nALL per-type rows reconcile to §3.1 full-panel n AND corrected African file: {all_ok}")
    return all_ok


if __name__ == "__main__":
    ok = main()
    (PROJ / "v5/00_logs/t1_origin_breakdown.log").write_text("\n".join(log) + "\n")
    print(f"\nLog -> {PROJ/'v5/00_logs/t1_origin_breakdown.log'}")
