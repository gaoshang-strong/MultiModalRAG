#!/usr/bin/env python3
"""
draw_pathway.py  —  Generalized pathway visualization for the neo4j-query skill.

Usage
-----
    python draw_pathway.py --config pathway.json --out figure.png [--dpi 150]

Config JSON format
------------------
{
  "title":    "ATM–TP53 Signal Flow",
  "subtitle": "optional note at bottom",   // optional
  "layout":   "flow",                       // "flow" | "columns"
  "figsize":  [15, 11],                     // optional, default auto
  "xlim":     [-8, 9.5],                    // optional, for flow only
  "ylim":     [-0.8, 10.2],                 // optional, for flow only

  // ── Zones (optional, for flow layout only) ──
  "zones": [
    {"comp": "left",   "label": "Left Pathway",     "x": [-8,-1],  "y": [-1,10]},
    {"comp": "right",  "label": "Right Pathway",    "x": [1.5,9],  "y": [3.5,10]},
    {"comp": "center", "label": "Convergence Zone", "x": [-2.5,4], "y": [-1,9.5], "linestyle": "--"}
  ],

  // ── Nodes ──
  // flow layout:    provide x, y coordinates (data units)
  // columns layout: provide "column": "left"|"center"|"right" instead of x/y
  "nodes": [
    {"key": "dna", "label": "DNA\nDamage", "comp": "trigger", "x": 0.0,  "y": 9.0},
    {"key": "atm", "label": "ATM",         "comp": "left",    "x": -3.5, "y": 6.0}
  ],

  // ── Edges ──
  // type: phos | tx | recruit | ppi | ubi_inh
  "edges": [
    {"src": "dna", "dst": "atm", "type": "phos", "label": "", "rad": 0.0}
  ]
}

Compartment color map (fixed)
------------------------------
  left    → blue   (#1B6BA8)
  right   → red    (#C0392B)
  center  → purple (#6A0DAD)
  trigger → dark   (#333333)
  output  → green  (#27AE60)

Edge type map (fixed)
---------------------
  phos     → orange solid arrow      (phosphorylation)
  tx       → green  dashed arrow     (transcription)
  recruit  → gray   dotted arrow     (recruitment / PPI)
  ppi      → light-gray solid arrow  (protein-protein interaction)
  ubi_inh  → red    solid blunt      (ubiquitin-mediated inhibition)
"""

import argparse
import json
import sys
import os
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ── Node visual constants ────────────────────────────────────────────────────
NW, NH = 1.45, 0.52    # node box width / height (data coords)
GAP    = 0.07           # gap between arrow tip and node edge

NFACE = {'left':'#CCDFF5','right':'#F9D0CA','center':'#E4D0F0',
         'trigger':'#E0E0E0','output':'#C8EDD0'}
NEDGE = {'left':'#1B6BA8','right':'#C0392B','center':'#6A0DAD',
         'trigger':'#333333','output':'#27AE60'}
BG_FACE = {'left':'#EBF5FF','right':'#FFF0EE','center':'#F5EEFF'}
BG_EDGE = {'left':'#2166AC','right':'#C0392B','center':'#6A0DAD'}
BG_LS   = {'left':'-','right':'-','center':'--'}

# ── Edge visual constants ────────────────────────────────────────────────────
ECOL = {'phos':'#E07B25','tx':'#2E8B57','recruit':'#7777AA',
        'ppi':'#AAAAAA','ubi_inh':'#C0392B'}
ELS  = {'phos':'-','tx':'--','recruit':':','ppi':'-','ubi_inh':'-'}
ELW  = {'phos':1.4,'tx':1.2,'recruit':1.0,'ppi':0.9,'ubi_inh':1.3}
EARR = {'phos':'->','tx':'->','recruit':'->','ppi':'->','ubi_inh':'-|>'}

LEGEND_LABELS = {
    'left':    'Left pathway',
    'right':   'Right pathway',
    'center':  'Convergence zone',
    'output':  'Output / cell fate',
}


# ── Drawing helpers ──────────────────────────────────────────────────────────

def _edge_pt(center, other):
    """Nearest point on node bounding box edge toward 'other'."""
    cx, cy = center
    dx, dy = other[0] - cx, other[1] - cy
    if abs(dx) < 1e-9 and abs(dy) < 1e-9:
        return center
    hw, hh = NW / 2 + GAP, NH / 2 + GAP
    t = min(hw / abs(dx) if abs(dx) > 1e-9 else 1e9,
            hh / abs(dy) if abs(dy) > 1e-9 else 1e9)
    return (cx + t * dx, cy + t * dy)


def draw_node(ax, label, pos, comp, fs=7.5):
    x, y = pos
    face = NFACE.get(comp, '#EEEEEE')
    edge = NEDGE.get(comp, '#444444')
    ax.add_patch(FancyBboxPatch(
        (x - NW / 2, y - NH / 2), NW, NH,
        boxstyle="round,pad=0.08",
        facecolor=face, edgecolor=edge,
        linewidth=1.8, zorder=3))
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fs, fontweight='bold', color=edge,
            multialignment='center', zorder=4)


def draw_edge(ax, p1, p2, etype, label='', rad=0.0, alpha=1.0):
    sp = _edge_pt(p1, p2)
    ep = _edge_pt(p2, p1)
    c  = ECOL.get(etype, '#888888')
    ax.annotate('', xy=ep, xytext=sp,
        arrowprops=dict(
            arrowstyle=EARR.get(etype, '->'),
            color=c, lw=ELW.get(etype, 1.0),
            linestyle=ELS.get(etype, '-'),
            connectionstyle=f'arc3,rad={rad}',
            alpha=alpha,
        ), zorder=2)
    if label:
        ax.text((sp[0] + ep[0]) / 2 + 0.1, (sp[1] + ep[1]) / 2,
                label, fontsize=5, color=c,
                ha='left', va='center', style='italic', zorder=5)


def draw_zones(ax, zones):
    """Draw background rectangles for named zones."""
    for z in zones:
        comp = z['comp']
        if comp not in BG_FACE:
            continue
        x0, x1 = z['x']
        y0, y1 = z['y']
        ls = z.get('linestyle', BG_LS.get(comp, '-'))
        ax.add_patch(plt.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            facecolor=BG_FACE[comp], edgecolor=BG_EDGE[comp],
            linewidth=1.0, alpha=0.28, zorder=0, linestyle=ls))
        ax.text(x0 + 0.15, y1 - 0.15, z.get('label', ''),
                fontsize=9, fontweight='bold',
                color=BG_EDGE[comp], va='top', alpha=0.75)


def draw_columns_zones(ax, xlim, ylim):
    """Draw 3 vertical bands for columns layout."""
    mid_x = 0.0
    x_split = (xlim[0] + mid_x) / 2 - 0.1  # ~-3
    x_split2 = (xlim[1] + mid_x) / 2 + 0.1  # ~+3
    zones = [
        {'comp': 'left',   'label': 'Left Pathway',     'x': [xlim[0], x_split],  'y': ylim},
        {'comp': 'center', 'label': 'Convergence Zone', 'x': [x_split, x_split2], 'y': ylim, 'linestyle': '--'},
        {'comp': 'right',  'label': 'Right Pathway',    'x': [x_split2, xlim[1]], 'y': ylim},
    ]
    draw_zones(ax, zones)


def make_legend(ax, bbox=(0.0, 0.0)):
    handles = [
        mpatches.Patch(facecolor=NFACE['left'],    edgecolor=NEDGE['left'],    label='Left pathway'),
        mpatches.Patch(facecolor=NFACE['right'],   edgecolor=NEDGE['right'],   label='Right pathway'),
        mpatches.Patch(facecolor=NFACE['center'],  edgecolor=NEDGE['center'],  label='Convergence zone'),
        mpatches.Patch(facecolor=NFACE['output'],  edgecolor=NEDGE['output'],  label='Output / cell fate'),
        plt.Line2D([0],[0], color=ECOL['phos'],    lw=1.4, label='Phosphorylation'),
        plt.Line2D([0],[0], color=ECOL['tx'],      lw=1.2, linestyle='--', label='Transcription'),
        plt.Line2D([0],[0], color=ECOL['recruit'], lw=1.0, linestyle=':',  label='Recruitment / PPI'),
        plt.Line2D([0],[0], color=ECOL['ubi_inh'], lw=1.3, label='Ubiquitin inhibition'),
    ]
    ax.legend(handles=handles, loc='lower left', fontsize=6.5, framealpha=0.9,
              bbox_to_anchor=bbox, ncol=2)


# ── Layout computation ───────────────────────────────────────────────────────

def compute_positions_flow(nodes):
    """Read x/y directly from each node definition."""
    pos = {}
    for n in nodes:
        if 'x' not in n or 'y' not in n:
            raise ValueError(f"Node '{n['key']}' missing x/y in flow layout")
        pos[n['key']] = (float(n['x']), float(n['y']))
    return pos


def compute_positions_columns(nodes,
                               x_left=-5.5, x_center=0.0, x_right=5.5,
                               y_top=9.0, y_bot=0.8,
                               y_trigger=10.2, y_output=-0.5):
    """
    Auto-distribute nodes in three columns by comp.
    Nodes with comp='trigger' are spread at the top; 'output' at the bottom.
    """
    buckets = {'left': [], 'center': [], 'right': [], 'trigger': [], 'output': []}
    for n in nodes:
        c = n.get('comp', 'center')
        buckets.setdefault(c, []).append(n['key'])

    pos = {}
    for comp, x in [('left', x_left), ('center', x_center), ('right', x_right)]:
        keys = buckets[comp]
        if not keys:
            continue
        ys = np.linspace(y_top, y_bot, len(keys)) if len(keys) > 1 else [( y_top + y_bot) / 2]
        for key, y in zip(keys, ys):
            pos[key] = (x, float(y))

    # Trigger nodes: spread horizontally at top
    trig = buckets['trigger']
    if trig:
        xs = np.linspace(x_left, x_right, len(trig)) if len(trig) > 1 else [0.0]
        for key, x in zip(trig, xs):
            pos[key] = (float(x), y_trigger)

    # Output nodes: spread horizontally at bottom
    out = buckets['output']
    if out:
        xs = np.linspace(x_left, x_right, len(out)) if len(out) > 1 else [0.0]
        for key, x in zip(out, xs):
            pos[key] = (float(x), y_output)

    return pos


# ── Auto axis limits from node positions ────────────────────────────────────

def auto_limits(pos, pad_x=1.5, pad_y=0.8):
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    return ([min(xs) - pad_x, max(xs) + pad_x],
            [min(ys) - pad_y, max(ys) + pad_y])


# ── Main render ──────────────────────────────────────────────────────────────

def render(config, outfile, dpi=150):
    layout   = config.get('layout', 'flow')
    nodes    = config['nodes']
    edges    = config.get('edges', [])
    title    = config.get('title', '')
    subtitle = config.get('subtitle', '')

    # Build position map
    if layout == 'columns':
        pos = compute_positions_columns(nodes)
    else:
        pos = compute_positions_flow(nodes)

    # Axis limits
    xlim_cfg = config.get('xlim')
    ylim_cfg = config.get('ylim')
    auto_xl, auto_yl = auto_limits(pos)
    xlim = xlim_cfg if xlim_cfg else auto_xl
    ylim = ylim_cfg if ylim_cfg else auto_yl

    # Figure size
    fw, fh = config.get('figsize', [15, 11])
    fig, ax = plt.subplots(figsize=(fw, fh))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis('off')

    # Background zones
    if layout == 'columns':
        draw_columns_zones(ax, xlim, ylim)
    elif 'zones' in config:
        draw_zones(ax, config['zones'])

    # Node lookup
    node_map = {n['key']: n for n in nodes}

    # Edges (draw first so nodes appear on top)
    edge_alpha = 0.70 if layout == 'columns' else 1.0
    for e in edges:
        src, dst = e['src'], e['dst']
        if src not in pos or dst not in pos:
            print(f"  WARNING: edge {src}→{dst} skipped (missing node)", file=sys.stderr)
            continue
        draw_edge(ax, pos[src], pos[dst],
                  etype=e.get('type', 'ppi'),
                  label=e.get('label', ''),
                  rad=float(e.get('rad', 0.0)),
                  alpha=edge_alpha)

    # Nodes
    fs = 6.8 if layout == 'columns' else 7.5
    for n in nodes:
        draw_node(ax, n['label'], pos[n['key']], n.get('comp', 'center'), fs=fs)

    # Legend
    make_legend(ax)

    # Title + subtitle
    if title:
        ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
    if subtitle:
        ax.text(0.5, -0.04, subtitle, transform=ax.transAxes,
                ha='center', va='top', fontsize=7.5,
                style='italic', color='#555555')

    fig.savefig(outfile, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outfile}')


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Draw a pathway diagram from a JSON config.')
    parser.add_argument('--config', required=True, help='Path to JSON config file')
    parser.add_argument('--out',    required=True, help='Output PNG path')
    parser.add_argument('--dpi',    type=int, default=150, help='Output DPI (default 150)')
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    render(config, args.out, dpi=args.dpi)


if __name__ == '__main__':
    main()
