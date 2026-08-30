#!/usr/bin/env python
"""
v5 METHODS_REPORT revision - TASK 2 render: annotated circular (fan) phylogeny of the
v5 panel (997 origin-normalized canonical genomes). Inner = ML phylogram (true branch
lengths, IQ-TREE2 --fast, GTR+F+G4). Outer rings annotate: HPV type (12 colours),
alpha-species clade (Alpha-9/7/6), and African origin (corrected geography, §3.1a/§7.4).
Per-type cluster labels, legend and scale bar; no per-tip text (997 tips).

Inputs: v5/03_phylogeny/v5_panel_iqtree.treefile, v5/03_phylogeny/v5_panel_metadata.tsv
Outputs: v5/06_reports/figures/phylo_tree_v5_panel.{png,svg,pdf}
"""
import math
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

P = Path(__file__).resolve().parents[1]  # repo root (design-scripts/ -> ..)
PHY = P / "v5/03_phylogeny"
FIGDIR = P / "v5/06_reports/figures"; FIGDIR.mkdir(parents=True, exist_ok=True)

# metadata
meta = {}
for ln in (PHY / "v5_panel_metadata.tsv").read_text().splitlines()[1:]:
    c = ln.split("\t")
    meta[c[0]] = dict(acc=c[1], type=c[2], clade=c[3], lineage=c[4], african=(c[5] == "Y"))

TYPES = ["HPV16","HPV31","HPV33","HPV35","HPV52","HPV58",   # Alpha-9
         "HPV18","HPV39","HPV45","HPV59",                    # Alpha-7
         "HPV56","HPV66"]                                    # Alpha-6
# 12 distinct type colours (grouped by clade hue family)
TYPE_COLOR = {
    "HPV16":"#08519c","HPV31":"#3182bd","HPV33":"#6baed6","HPV35":"#9ecae1","HPV52":"#2171b5","HPV58":"#c6dbef",
    "HPV18":"#a50f15","HPV39":"#de2d26","HPV45":"#fb6a4a","HPV59":"#fcae91",
    "HPV56":"#238b45","HPV66":"#74c476",
}
CLADE_COLOR = {"Alpha9":"#1f3b73","Alpha7":"#99000d","Alpha6":"#00441b"}
AFR_COLOR = "#e6a000"

tree = Phylo.read(str(PHY / "v5_panel_iqtree.treefile"), "newick")
tree.ladderize()
leaves = tree.get_terminals()
N = len(leaves)
parents = {}
for cl in tree.find_clades():
    for ch in cl.clades:
        parents[ch] = cl
ang = {}
for i, lf in enumerate(leaves):
    ang[lf] = 2 * math.pi * i / N
for cl in tree.find_clades(order="postorder"):
    if not cl.is_terminal():
        ang[cl] = sum(ang[c] for c in cl.clades) / len(cl.clades)
depth = tree.depths()           # branch-length distance from root
if max(depth.values()) == 0:    # safety: fall back to topology depths
    depth = tree.depths(unit_branch_lengths=True)
Rmax = max(depth[lf] for lf in leaves)

def xy(r, a):
    return r * math.cos(a), r * math.sin(a)

fig, ax = plt.subplots(figsize=(17, 17))
ax.set_aspect("equal"); ax.axis("off")

# --- branches (gray phylogram) ---
for cl in tree.find_clades():
    if cl in parents:
        p = parents[cl]; a = ang[cl]
        x0, y0 = xy(depth[p], a); x1, y1 = xy(depth[cl], a)
        ax.plot([x0, x1], [y0, y1], color="#555555", lw=0.35, zorder=1)
    if not cl.is_terminal():
        a0 = min(ang[c] for c in cl.clades); a1 = max(ang[c] for c in cl.clades)
        r = depth[cl]
        aa = [a0 + (a1 - a0) * k / 40 for k in range(41)]
        ax.plot([r * math.cos(t) for t in aa], [r * math.sin(t) for t in aa],
                color="#555555", lw=0.35, zorder=1)

# --- annotation rings ---
r_type = Rmax * 1.06
r_clade = Rmax * 1.11
r_afr = Rmax * 1.16

# faint OUTER radial guide (only the annulus just inside the ring) so shallow tips read
# cleanly against their ring tick without flooding the centre
guide_in = Rmax * 0.97
for lf in leaves:
    a = ang[lf]
    x0, y0 = xy(guide_in, a); x1, y1 = xy(r_type * 0.995, a)
    ax.plot([x0, x1], [y0, y1], color="#ededed", lw=0.15, zorder=0)
dw = math.pi / N * 1.15          # angular half-width of a tick (slight overlap -> continuous band)
for lf in leaves:
    m = meta.get(lf.name)
    if not m:
        continue
    a = ang[lf]
    # type ring
    aa = [a - dw, a + dw]
    ax.plot([r_type * math.cos(t) for t in aa], [r_type * math.sin(t) for t in aa],
            color=TYPE_COLOR.get(m["type"], "#cccccc"), lw=3.2, solid_capstyle="butt", zorder=2)
    # clade ring
    ax.plot([r_clade * math.cos(t) for t in aa], [r_clade * math.sin(t) for t in aa],
            color=CLADE_COLOR.get(m["clade"], "#cccccc"), lw=3.2, solid_capstyle="butt", zorder=2)
    # african marker
    if m["african"]:
        x, y = xy(r_afr, a)
        ax.plot([x], [y], marker="o", ms=2.4, color=AFR_COLOR, mec="none", zorder=3)

# --- per-type cluster labels at median angle ---
by_type = defaultdict(list)
for lf in leaves:
    m = meta.get(lf.name)
    if m: by_type[m["type"]].append(ang[lf])
r_lab = Rmax * 1.22
for ty, angs in by_type.items():
    a = sorted(angs)[len(angs) // 2]
    x, y = xy(r_lab, a)
    deg = math.degrees(a)
    rot = deg - 90 if -90 <= ((deg + 90) % 360) - 180 < 90 else deg + 90
    ha = "left" if math.cos(a) >= 0 else "right"
    ax.text(x, y, ty, fontsize=22, color=TYPE_COLOR[ty], fontweight="bold",
            ha=ha, va="center", rotation=deg if math.cos(a) >= 0 else deg + 180,
            rotation_mode="anchor", zorder=4)

# --- scale bar (substitutions/site) ---
sb = 0.05
x0 = -Rmax * 1.15; y0 = -Rmax * 1.30
ax.plot([x0, x0 + sb], [y0, y0], color="black", lw=2)
ax.text(x0 + sb / 2, y0 - Rmax * 0.03, f"{sb} subs/site", ha="center", va="top", fontsize=18)

# --- legends (figure-level, guaranteed visible to the right of the plot) ---
type_handles = [Line2D([0], [0], color=TYPE_COLOR[t], lw=6) for t in TYPES]
leg1 = fig.legend(type_handles, TYPES, title="HPV type (inner ring)",
                  loc="upper left", bbox_to_anchor=(0.64, 0.97),
                  fontsize=15, title_fontsize=16, frameon=True, facecolor="white", edgecolor="#cccccc")
fig.add_artist(leg1)
clade_handles = [Line2D([0], [0], color=CLADE_COLOR[c], lw=6) for c in ["Alpha9","Alpha7","Alpha6"]]
afr_handle = [Line2D([0], [0], marker="o", color="w", markerfacecolor=AFR_COLOR, ms=9)]
leg2 = fig.legend(clade_handles + afr_handle,
                  ["Alpha-9","Alpha-7","Alpha-6","African origin (corrected)"],
                  title="Clade (middle ring) / origin (outer dots)",
                  loc="upper left", bbox_to_anchor=(0.64, 0.52),
                  fontsize=15, title_fontsize=16, frameon=True, facecolor="white", edgecolor="#cccccc")
fig.add_artist(leg2)

ax.set_title("ML phylogeny of the v5 panel (997 canonical genomes) — "
             "IQ-TREE2 2.3.6 --fast, GTR+F+G4\n"
             "inner ring = HPV type; middle ring = alpha-clade; outer dots = African origin",
             fontsize=16)
lim = Rmax * 1.50
ax.set_xlim(-lim, lim); ax.set_ylim(-lim * 1.05, lim)
fig.subplots_adjust(left=0.01, right=0.62, top=0.93, bottom=0.04)

for ext in ("png", "svg", "pdf"):
    fig.savefig(FIGDIR / f"phylo_tree_v5_panel.{ext}",
                dpi=160 if ext == "png" else None)
print("wrote", [str(FIGDIR / f"phylo_tree_v5_panel.{e}") for e in ("png","svg","pdf")])
print(f"tips={N} Rmax={Rmax:.3f} African_tips={sum(1 for lf in leaves if meta.get(lf.name,{}).get('african'))}")
