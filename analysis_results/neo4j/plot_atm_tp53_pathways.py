#!/usr/bin/env python3
"""
ATM-TP53 Pathway Network Visualization — standalone exploratory figures
Produces:
  20260606_ATM_TP53_flow.png     — signal flow diagram (~23 nodes)
  20260606_ATM_TP53_network.png  — three-compartment detailed network (~43 nodes)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os

OUT = os.path.dirname(os.path.abspath(__file__))

# ── Aesthetics ──────────────────────────────────────────────────────────────
NW, NH = 1.45, 0.52
GAP    = 0.07

NFACE = {'atm':'#CCDFF5','tp53':'#F9D0CA','conv':'#E4D0F0',
         'dna':'#E0E0E0', 'outcome':'#C8EDD0'}
NEDGE = {'atm':'#1B6BA8','tp53':'#C0392B','conv':'#6A0DAD',
         'dna':'#333333', 'outcome':'#27AE60'}
ECOL  = {'phos':'#E07B25','tx':'#2E8B57','recruit':'#7777AA',
         'ppi':'#AAAAAA', 'ubi_inh':'#C0392B'}
ELS   = {'phos':'-','tx':'--','recruit':':','ppi':'-','ubi_inh':'-'}
ELW   = {'phos':1.4,'tx':1.2,'recruit':1.0,'ppi':0.9,'ubi_inh':1.3}
EARR  = {'phos':'->','tx':'->','recruit':'->','ppi':'->','ubi_inh':'-|>'}


def edge_pt(center, other, hw=None, hh=None):
    hw = hw or NW/2 + GAP
    hh = hh or NH/2 + GAP
    cx, cy = center
    dx, dy = other[0]-cx, other[1]-cy
    if abs(dx) < 1e-9 and abs(dy) < 1e-9:
        return center
    t = min(hw/abs(dx) if abs(dx) > 1e-9 else 1e9,
            hh/abs(dy) if abs(dy) > 1e-9 else 1e9)
    return (cx + t*dx, cy + t*dy)


def draw_node(ax, label, pos, comp, fs=7.5):
    x, y = pos
    ax.add_patch(FancyBboxPatch(
        (x-NW/2, y-NH/2), NW, NH,
        boxstyle="round,pad=0.08",
        facecolor=NFACE[comp], edgecolor=NEDGE[comp],
        linewidth=1.8, zorder=3))
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fs, fontweight='bold',
            color=NEDGE[comp], multialignment='center', zorder=4)


def draw_edge(ax, p1, p2, etype, label='', rad=0.0, alpha=1.0, hw=None, hh=None):
    sp = edge_pt(p1, p2, hw, hh)
    ep = edge_pt(p2, p1, hw, hh)
    c  = ECOL[etype]
    ax.annotate('', xy=ep, xytext=sp,
        arrowprops=dict(
            arrowstyle=EARR[etype], color=c,
            lw=ELW[etype], linestyle=ELS[etype],
            connectionstyle=f'arc3,rad={rad}',
            alpha=alpha,
        ), zorder=2)
    if label:
        ax.text((sp[0]+ep[0])/2+0.1, (sp[1]+ep[1])/2, label,
                fontsize=5, color=c, ha='left', va='center',
                style='italic', zorder=5)


def make_legend(ax, loc='lower left', bbox=(0.0, 0.0)):
    handles = [
        mpatches.Patch(facecolor=NFACE['atm'],     edgecolor=NEDGE['atm'],     label='ATM-specific'),
        mpatches.Patch(facecolor=NFACE['tp53'],    edgecolor=NEDGE['tp53'],    label='TP53-specific'),
        mpatches.Patch(facecolor=NFACE['conv'],    edgecolor=NEDGE['conv'],    label='Convergence zone'),
        mpatches.Patch(facecolor=NFACE['outcome'], edgecolor=NEDGE['outcome'], label='Cell fate output'),
        plt.Line2D([0],[0], color=ECOL['phos'],    lw=1.4, label='Phosphorylation'),
        plt.Line2D([0],[0], color=ECOL['tx'],      lw=1.2, linestyle='--', label='Transcription'),
        plt.Line2D([0],[0], color=ECOL['recruit'], lw=1.0, linestyle=':', label='Recruitment / PPI'),
        plt.Line2D([0],[0], color=ECOL['ubi_inh'], lw=1.3, label='Ubiquitin inhibition'),
    ]
    ax.legend(handles=handles, loc=loc, fontsize=6.5, framealpha=0.9,
              bbox_to_anchor=bbox, ncol=2)


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Signal flow diagram (hierarchical)
# ════════════════════════════════════════════════════════════════════════════

F1_NODES = {
    'dna_dsb':    ('DNA\nDamage',          'dna',      (0.0,  9.0)),
    'mrn':        ('MRN\nComplex',         'atm',      (-3.5, 7.5)),
    'kat5':       ('KAT5\n(TIP60)',        'conv',     (0.8,  7.5)),
    'atm':        ('ATM',                  'atm',      (-3.5, 6.0)),
    'tp53':       ('TP53',                 'tp53',     (3.5,  6.0)),
    'mdm2':       ('MDM2',                 'tp53',     (6.5,  6.0)),
    'h2ax':       ('γH2AX',                'atm',      (-6.0, 4.5)),
    'mdc1':       ('MDC1',                 'atm',      (-4.0, 4.5)),
    'chek2':      ('CHEK2',                'conv',     (-1.5, 4.5)),
    'atr':        ('ATR',                  'conv',     (0.5,  4.5)),
    'cdkn1a':     ('CDKN1A\n(p21)',        'tp53',     (2.5,  4.5)),
    'bax_puma':   ('BAX/PUMA',             'tp53',     (5.0,  4.5)),
    'gadd45':     ('GADD45',               'tp53',     (7.5,  4.5)),
    'rnf8':       ('RNF8',                 'atm',      (-5.5, 3.0)),
    'chek1':      ('CHEK1',                'conv',     (-0.5, 3.0)),
    'brca1':      ('BRCA1',                'conv',     (2.5,  3.0)),
    'rnf168':     ('RNF168',               'atm',      (-5.5, 1.5)),
    'bp53bp1':    ('53BP1',                'atm',      (-3.5, 1.5)),
    'rad51':      ('RAD51',                'conv',     (2.5,  1.5)),
    'arrest':     ('Cell Cycle\nArrest',   'outcome',  (-4.0, 0.0)),
    'repair':     ('DNA\nRepair',          'outcome',  (0.0,  0.0)),
    'apoptosis':  ('Apoptosis',            'outcome',  (3.5,  0.0)),
    'senescence': ('Senescence',           'outcome',  (7.0,  0.0)),
}

F1_EDGES = [
    # src, dst, type, label, rad
    ('dna_dsb',  'mrn',       'ppi',     '',              0.0),
    ('dna_dsb',  'kat5',      'ppi',     '',              0.2),
    ('mrn',      'atm',       'recruit', 'activates ATM', 0.0),
    ('kat5',     'atm',       'phos',    'K3016 acetyl',  0.2),
    ('kat5',     'tp53',      'phos',    'K382 acetyl',   0.25),
    ('atm',      'tp53',      'phos',    'S15 stabilize', 0.0),
    ('atm',      'mdm2',      'phos',    'S395 inhibit',  0.2),
    ('mdm2',     'tp53',      'ubi_inh', 'degrades',      0.0),
    ('atm',      'h2ax',      'phos',    'γH2AX mark',    0.0),
    ('atm',      'mdc1',      'phos',    '',              0.0),
    ('atm',      'chek2',     'phos',    '',              0.0),
    ('atm',      'atr',       'ppi',     '',              0.0),
    ('h2ax',     'mdc1',      'ppi',     '',              0.0),
    ('mdc1',     'rnf8',      'recruit', '',              0.0),
    ('rnf8',     'rnf168',    'phos',    'H2A-Ub cascade',0.0),
    ('rnf168',   'bp53bp1',   'recruit', '',              0.0),
    ('atr',      'chek1',     'phos',    '',              0.0),
    ('chek2',    'tp53',      'phos',    '',              0.3),
    ('bp53bp1',  'chek2',     'ppi',     '',              0.25),
    ('tp53',     'cdkn1a',    'tx',      '',              0.0),
    ('tp53',     'bax_puma',  'tx',      '',              0.0),
    ('tp53',     'gadd45',    'tx',      '',              0.0),
    ('tp53',     'brca1',     'tx',      '',              0.35),
    ('atm',      'brca1',     'phos',    '',              0.15),
    ('brca1',    'rad51',     'recruit', 'HR repair',     0.0),
    # Outcomes
    ('chek1',    'arrest',    'phos',    '',              0.0),
    ('chek2',    'arrest',    'phos',    '',              0.2),
    ('cdkn1a',   'arrest',    'ppi',     '',              0.0),
    ('rad51',    'repair',    'ppi',     '',              0.0),
    ('brca1',    'repair',    'ppi',     '',              0.2),
    ('bp53bp1',  'repair',    'ppi',     '',              0.2),
    ('bax_puma', 'apoptosis', 'tx',      '',              0.0),
    ('chek2',    'apoptosis', 'phos',    '',              0.35),
    ('gadd45',   'senescence','tx',      '',              0.0),
]


def fig1_flow():
    fig, ax = plt.subplots(figsize=(15, 11))
    ax.set_xlim(-8.0, 9.5)
    ax.set_ylim(-0.8, 10.2)
    ax.axis('off')

    # Background zones
    ax.add_patch(plt.Rectangle((-8.0, -0.6), 6.8, 10.5,
        facecolor='#EBF5FF', edgecolor='#2166AC', lw=1.0,
        alpha=0.30, zorder=0))
    ax.add_patch(plt.Rectangle((1.5, 3.5), 8.3, 6.8,
        facecolor='#FFF0EE', edgecolor='#C0392B', lw=1.0,
        alpha=0.30, zorder=0))
    ax.add_patch(plt.Rectangle((-2.5, -0.6), 6.5, 9.5,
        facecolor='#F5EEFF', edgecolor='#6A0DAD', lw=1.0,
        alpha=0.20, zorder=0, linestyle='--'))

    ax.text(-7.6, 9.8, 'ATM Kinase Cascade', fontsize=9, fontweight='bold',
            color='#2166AC', va='top', alpha=0.75)
    ax.text(9.2, 10.1, 'TP53 Transcriptional\nProgramme', fontsize=9,
            fontweight='bold', color='#C0392B', va='top', ha='right', alpha=0.75)
    ax.text(-2.2, 9.7, 'Convergence Zone', fontsize=8, fontweight='bold',
            color='#6A0DAD', va='top', ha='left', alpha=0.65, style='italic')

    pos = {k: v[2] for k, v in F1_NODES.items()}
    for src, dst, etype, lbl, rad in F1_EDGES:
        draw_edge(ax, pos[src], pos[dst], etype, lbl, rad)
    for key, (label, comp, p) in F1_NODES.items():
        draw_node(ax, label, p, comp)

    make_legend(ax, loc='lower left', bbox=(0.0, 0.0))

    ax.set_title(
        'ATM–TP53 Signal Flow  ·  Mutual Exclusivity Mechanism in Gastric Cancer',
        fontsize=12, fontweight='bold', pad=8)

    note = ('Mutual exclusivity arises because ATM kinase cascade (blue) and TP53 '
            'transcriptional programme (red)\nboth converge on shared effectors '
            '(purple). Single-pathway loss is tolerated; dual loss silences all DDR output.')
    ax.text(0.5, -0.06, note, transform=ax.transAxes,
            ha='center', va='top', fontsize=7.5, style='italic', color='#555555')

    out = os.path.join(OUT, '20260606_ATM_TP53_flow.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out}')


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Three-compartment detailed network (column layout)
# ════════════════════════════════════════════════════════════════════════════

_atm_col  = ['ATM', 'MRE11', 'NBN', 'RAD50', 'MDC1', 'H2AX',
             'RNF8', 'RNF168', '53BP1', 'ABRAXAS1', 'UIMC1',
             'BABAM1', 'HERC2', 'BRCC3']
_conv_col = ['KAT5', 'CHEK2', 'CHEK1', 'ATR', 'BRCA1', 'RAD51',
             'PCNA', 'FOXO3', 'DYRK2', 'MAPK8', 'SETD2',
             'STK11', 'BLM', 'WRN', 'MSH2']
_tp53_col = ['TP53', 'MDM2', 'MDM4', 'CDKN1A', 'BAX',
             'PUMA', 'GADD45', 'CCNG1', 'XPC', 'E2F1']
_out_col  = ['Cell Cycle\nArrest', 'DNA\nRepair', 'Apoptosis', 'Senescence']

def _col_pos(names, x, y_top=9.0, y_bot=0.8):
    ys = np.linspace(y_top, y_bot, len(names))
    return {n: (x, float(y)) for n, y in zip(names, ys)}

F2_POS = {}
F2_POS.update(_col_pos(_atm_col,  -5.5))
F2_POS.update(_col_pos(_conv_col,  0.0))
F2_POS.update(_col_pos(_tp53_col,  5.5))
F2_POS.update({
    'Cell Cycle\nArrest': (-4.5, -0.5),
    'DNA\nRepair':         (-1.0, -0.5),
    'Apoptosis':           ( 2.5, -0.5),
    'Senescence':          ( 6.0, -0.5),
})

F2_COMP = {}
for n in _atm_col:  F2_COMP[n] = 'atm'
for n in _conv_col: F2_COMP[n] = 'conv'
for n in _tp53_col: F2_COMP[n] = 'tp53'
for n in _out_col:  F2_COMP[n] = 'outcome'

F2_EDGES = [
    # Within ATM column
    ('ATM',      'MRE11',    'ppi',     0.0),
    ('ATM',      'NBN',      'ppi',     0.0),
    ('ATM',      'RAD50',    'ppi',     0.0),
    ('ATM',      'MDC1',     'phos',    0.0),
    ('ATM',      'H2AX',     'phos',    0.0),
    ('H2AX',     'MDC1',     'ppi',     0.2),
    ('MDC1',     'RNF8',     'recruit', 0.0),
    ('RNF8',     'RNF168',   'phos',    0.0),
    ('RNF168',   '53BP1',    'recruit', 0.0),
    ('HERC2',    'RNF168',   'ppi',     0.25),
    ('53BP1',    'ABRAXAS1', 'ppi',     0.0),
    ('ABRAXAS1', 'UIMC1',    'ppi',     0.0),
    ('UIMC1',    'BABAM1',   'ppi',     0.0),
    # ATM → Convergence
    ('ATM',      'CHEK2',    'phos',    0.0),
    ('ATM',      'ATR',      'ppi',     0.0),
    ('ATM',      'BRCA1',    'phos',    0.0),
    ('ATM',      'KAT5',     'ppi',     0.2),
    ('MRE11',    'CHEK2',    'ppi',     0.2),
    ('53BP1',    'CHEK2',    'ppi',     0.0),
    # ATM → TP53
    ('ATM',      'TP53',     'phos',    0.0),
    # Convergence internal
    ('KAT5',     'CHEK2',    'ppi',     0.2),
    ('ATR',      'CHEK1',    'phos',    0.0),
    ('CHEK2',    'CHEK1',    'ppi',     0.2),
    ('BRCA1',    'RAD51',    'recruit', 0.0),
    ('BLM',      'RAD51',    'ppi',     0.0),
    ('WRN',      'RAD51',    'ppi',     0.0),
    ('PCNA',     'RAD51',    'ppi',     0.2),
    ('MSH2',     'ATR',      'recruit', 0.2),
    ('STK11',    'CHEK2',    'phos',    0.2),
    # Convergence → TP53
    ('KAT5',     'TP53',     'phos',    0.0),
    ('CHEK2',    'TP53',     'phos',    0.0),
    ('DYRK2',    'TP53',     'phos',    0.0),
    ('MAPK8',    'TP53',     'phos',    0.0),
    ('FOXO3',    'CDKN1A',   'tx',      0.0),
    # Within TP53 column
    ('MDM2',     'TP53',     'ubi_inh', 0.2),
    ('MDM4',     'TP53',     'ubi_inh', 0.2),
    ('E2F1',     'TP53',     'ppi',     0.2),
    ('TP53',     'MDM2',     'tx',      0.3),
    ('TP53',     'CDKN1A',   'tx',      0.0),
    ('TP53',     'BAX',      'tx',      0.0),
    ('TP53',     'PUMA',     'tx',      0.0),
    ('TP53',     'GADD45',   'tx',      0.0),
    ('TP53',     'CCNG1',    'tx',      0.0),
    ('TP53',     'XPC',      'tx',      0.0),
    # → Outcomes
    ('CHEK1',    'Cell Cycle\nArrest',  'phos', 0.0),
    ('CHEK2',    'Cell Cycle\nArrest',  'phos', 0.2),
    ('CDKN1A',   'Cell Cycle\nArrest',  'ppi',  0.0),
    ('RAD51',    'DNA\nRepair',         'ppi',  0.0),
    ('BRCA1',    'DNA\nRepair',         'ppi',  0.2),
    ('53BP1',    'DNA\nRepair',         'ppi',  0.2),
    ('XPC',      'DNA\nRepair',         'ppi',  0.0),
    ('BAX',      'Apoptosis',           'tx',   0.0),
    ('PUMA',     'Apoptosis',           'tx',   0.0),
    ('CHEK2',    'Apoptosis',           'phos', 0.3),
    ('GADD45',   'Senescence',          'tx',   0.0),
    ('CCNG1',    'Senescence',          'tx',   0.2),
]


def fig2_network():
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_xlim(-7.8, 8.5)
    ax.set_ylim(-1.2, 10.5)
    ax.axis('off')

    # Column background bands
    ax.add_patch(plt.Rectangle((-7.8, -1.0), 4.8, 11.2,
        facecolor='#EBF5FF', edgecolor='#2166AC', lw=1.0, alpha=0.25, zorder=0))
    ax.add_patch(plt.Rectangle((-2.8, -1.0), 5.6, 11.2,
        facecolor='#F5EEFF', edgecolor='#6A0DAD', lw=1.0, alpha=0.20, zorder=0, linestyle='--'))
    ax.add_patch(plt.Rectangle((3.0, -1.0), 5.5, 11.2,
        facecolor='#FFF0EE', edgecolor='#C0392B', lw=1.0, alpha=0.25, zorder=0))

    ax.text(-7.5, 10.2, 'ATM Cascade\n(DSB repair / checkpoint)',
            fontsize=9, fontweight='bold', color='#1B6BA8', va='top')
    ax.text( 0.0, 10.2, 'Convergence Zone\n(shared effectors)',
            fontsize=9, fontweight='bold', color='#6A0DAD', va='top', ha='center')
    ax.text( 8.2, 10.2, 'TP53 Programme\n(transcription / apoptosis)',
            fontsize=9, fontweight='bold', color='#C0392B', va='top', ha='right')

    # Draw edges
    for rec in F2_EDGES:
        src, dst, etype, rad = rec
        draw_edge(ax, F2_POS[src], F2_POS[dst], etype, rad=rad, alpha=0.65,
                  hw=NW/2+GAP+0.05, hh=NH/2+GAP+0.02)

    # Draw nodes (small font for detail figure)
    for name, comp in F2_COMP.items():
        draw_node(ax, name, F2_POS[name], comp, fs=6.8)

    make_legend(ax, loc='lower left', bbox=(0.0, 0.0))

    ax.set_title(
        'ATM–TP53 Detailed Network  ·  Three-Compartment Model (~43 nodes)',
        fontsize=12, fontweight='bold', pad=8)

    out = os.path.join(OUT, '20260606_ATM_TP53_network.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out}')


if __name__ == '__main__':
    fig1_flow()
    fig2_network()
    print('Done.')
