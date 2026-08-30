#!/usr/bin/env python
"""
METHODS_REPORT addition - tiling/coverage plot of the FINAL v5 scheme.
Visualisation only (reads existing design coordinates; no re-design/analysis).

For each of the 12 types, draws the genome as a horizontal axis (0..L canonical bp)
and the 6 amplicons as ARTIC-style staggered tiles: Pool 1 above the genome line,
Pool 2 below. The origin-junction amplicon (left_canstart > right_canend) is drawn
wrapping the origin as two segments. Overlaps between adjacent (alternate-pool)
amplicons are visible as horizontal extent shared across the two rows.

Inputs : v5/04_design/HPV{t}.amplicons.tsv, v5/06_reports/coverage_map.tsv
Outputs: v5/06_reports/figures/amplicon_tiling_v5.{png,svg,pdf}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

P = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
DES = P / "v5/04_design"
FIG = P / "v5/06_reports/figures"; FIG.mkdir(parents=True, exist_ok=True)
TYPES = [16, 18, 31, 33, 35, 39, 45, 52, 56, 58, 59, 66]
POOL_COLOR = {1: "#2c7fb8", 2: "#e6802e"}

# genome lengths + coverage from coverage_map.tsv
glen = {}; cov = {}
for ln in (P / "v5/06_reports/coverage_map.tsv").read_text().splitlines()[1:]:
    c = ln.split("\t"); glen[c[0]] = int(c[1]); cov[c[0]] = c[3]

fig, axes = plt.subplots(12, 1, figsize=(12, 18))
for ax, t in zip(axes, TYPES):
    ty = f"HPV{t}"; L = glen[ty]
    amps = []
    for ln in (DES / f"{ty}.amplicons.tsv").read_text().splitlines()[1:]:
        c = ln.split("\t")
        amps.append(dict(amp=c[0], pool=int(c[1]), prod=int(c[2]),
                         s=int(c[3]), e=int(c[4]), junc=(c[5] == "Y")))
    # genome backbone
    ax.hlines(0, 0, L, color="#444444", lw=1.4, zorder=1)
    ax.plot([0], [0], marker="|", ms=14, color="black", zorder=3)   # origin tick
    h = 0.7
    for a in amps:
        y = 0.55 if a["pool"] == 1 else -0.55 - h
        col = POOL_COLOR[a["pool"]]
        ec = "black" if a["junc"] else "none"
        lw = 1.3 if a["junc"] else 0
        hatch = "////" if a["junc"] else None
        segs = [(a["s"], a["e"])] if a["s"] <= a["e"] else [(a["s"], L), (0, a["e"])]
        for (xs, xe) in segs:
            ax.add_patch(Rectangle((xs, y), xe - xs, h, facecolor=col, edgecolor=ec,
                                   lw=lw, hatch=hatch, alpha=0.85, zorder=2))
        # amplicon label at midpoint of the (first) segment
        mx = (segs[0][0] + segs[0][1]) / 2
        ax.text(mx, y + h / 2, a["amp"], ha="center", va="center", fontsize=7,
                color="white", fontweight="bold", zorder=4)
    ax.set_xlim(-L * 0.01, L * 1.01)
    ax.set_ylim(-1.7, 1.5)
    ax.set_yticks([])
    ax.set_ylabel(ty, rotation=0, ha="right", va="center", fontsize=10, fontweight="bold")
    ax.set_title(f"{ty}  —  {L:,} bp  ·  {len(amps)} amplicons  ·  {cov[ty]}% coverage incl. origin",
                 fontsize=9, loc="left")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    if ax is not axes[-1]:
        ax.set_xlabel("")
    else:
        ax.set_xlabel("canonical genome position (bp)", fontsize=9)

# legend
handles = [Rectangle((0, 0), 1, 1, facecolor=POOL_COLOR[1]),
           Rectangle((0, 0), 1, 1, facecolor=POOL_COLOR[2]),
           Rectangle((0, 0), 1, 1, facecolor="#bbbbbb", edgecolor="black", hatch="////"),
           Line2D([0], [0], marker="|", ms=12, color="black", lw=0)]
labels = ["Pool 1 (above line)", "Pool 2 (below line)",
          "origin-junction amplicon (wraps position 0)", "circular origin (position 0)"]
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.005),
           ncol=4, fontsize=9, frameon=False)
fig.suptitle("Final v5 tiling scheme — 6 amplicons / 2 alternating pools per type, "
             "100% coverage including the circular origin", fontsize=13, y=1.02)
fig.tight_layout(rect=[0, 0, 1, 0.99])
for ext in ("png", "svg", "pdf"):
    fig.savefig(FIG / f"amplicon_tiling_v5.{ext}", bbox_inches="tight",
                dpi=160 if ext == "png" else None)
print("wrote", [str(FIG / f"amplicon_tiling_v5.{e}") for e in ("png", "svg", "pdf")])
