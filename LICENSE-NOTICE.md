# Licensing notice — dual license by file type

This repository is **dual-licensed**, split by file type:

- **Code → MIT** — see [`LICENSE-MIT.txt`](LICENSE-MIT.txt).
  Covers all scripts and source: `*.py`, `*.sh`, `Snakefile`, `*.nf`, `*.lua`,
  `*.yml`/`*.yaml`, `nextflow.config`, and similar executable/source files.
- **Data & documents → CC-BY-4.0** — see [`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt).
  Covers the primer-scheme outputs and design data (`primer.bed`, `reference.fasta`, the
  `*.tsv` tables — e.g. `primer_order_sheet_v5.tsv`, `scheme_summary.tsv`,
  `coordinate_rederivation.tsv`, and the inclusivity / specificity / lineage / pool tables —
  plus `MANIFEST.txt`, `SHA256SUMS`, `index.json`), the reports (`METHODS_REPORT.*`,
  `FINAL_REPORT_*`), and all figures.

**CC-BY-4.0 attribution.** Bright Adu, ORCID
[0000-0002-3622-2137](https://orcid.org/0000-0002-3622-2137), Genomics and Bioinformatics
Core, Noguchi Memorial Institute for Medical Research, University of Ghana, Legon, Ghana.
Please **cite via [`CITATION.cff`](CITATION.cff)**.

## Genuinely ambiguous files (assignment stated)
- `*.yml` conda env specs, `nextflow.config`, `*.lua`, and `preamble.tex` — build / source
  config; licensed as **code (MIT)**.
- `sample_sheet_TEMPLATE.csv` — a data template (not real data); licensed **CC-BY-4.0**.
- `index.json`, `MANIFEST.txt`, `SHA256SUMS`, `*.fai` — data / metadata describing the
  scheme; licensed **CC-BY-4.0**.
- `README.md` and this `LICENSE-NOTICE.md` are documentation (**CC-BY-4.0**); `CITATION.cff`
  and `.gitignore` are repository metadata.

> Note: GitHub may show **"View license"** rather than a single SPDX badge for a
> dual-licensed repository — this is expected, not an error.
