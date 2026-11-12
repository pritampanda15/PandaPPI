# pandappi/visualization/interaction_map.py
"""
Publication figures for PandaPPI: interface maps, interface reports, atom-level
contact views and the summary dashboard. All are drawn from the tidy interaction
table made by reports.generator.create_interaction_report().

Style follows Nature figure conventions: 183 mm wide, Helvetica/Arial, 5-7 pt
text, 8 pt bold lowercase panel letters, NPG palette, PNG (300 dpi) plus vector PDF.
"""

import functools
import os
import re
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

WIDTH = 7.2  # inches = 183 mm, a Nature double column
INK = '#424242'
TEAL, RED, NAVY, BLUE, ORANGE, DARK, GREY = '#00A087', '#E64B35', '#3C5488', '#4DBBD5', '#F39B7F', '#3D1B24', '#B0B0B0'

NATURE_RC = {
    # Arial before Helvetica: matplotlib cannot reach the bold face inside macOS's Helvetica.ttc
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 6, 'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 5.5, 'ytick.labelsize': 6,
    'legend.fontsize': 5.5, 'legend.title_fontsize': 6, 'legend.frameon': False,
    'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2, 'ytick.major.size': 2, 'lines.linewidth': 1.0, 'grid.linewidth': 0.3,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'figure.facecolor': 'white',
}


def nature(func):
    """Draw with the Nature style without touching the caller's matplotlib settings."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with plt.rc_context(NATURE_RC):
            return func(*args, **kwargs)
    return wrapper


# label -> (colour, line style, badge). Colour = chemistry (teal polar, red charged, navy pi,
# charcoal apolar, dark covalent/metal); line style and badge tell the types of one family apart.
TYPE_STYLE = {
    'Hydrogen bond': (TEAL, '-', 'H'),
    'Salt bridge / ionic': (RED, '-', 'I'),
    'Disulfide': (DARK, '-', 'SS'),
    'Pi-pi stacking': (NAVY, '-', 'ππ'),
    'Cation-pi': (RED, '--', '+π'),
    'CH-pi': (NAVY, '--', 'Cπ'),
    'Hydrophobic': (INK, ':', 'h'),
    'Water-mediated': (TEAL, '--', 'W'),
    'Metal coordination': (DARK, '--', 'M'),
    'Halogen bond': (ORANGE, '-', 'X'),
    'Sulfur-aromatic': (NAVY, '-.', 'Sπ'),
    'Amide-aromatic': (NAVY, ':', 'Aπ'),
    'Amide-amide': (TEAL, ':', 'A'),
    'Van der Waals': (GREY, ':', 'v'),
}
PEPTIDE_MAX = 30   # chains up to this many residues are drawn as a peptide (radial map)
BURIED_MIN = 40.0  # Å² a residue must bury to get the heavy outline

RESIDUE_CLASS = {
    **dict.fromkeys(['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'GLY'], 'Hydrophobic'),
    **dict.fromkeys(['PHE', 'TYR', 'TRP'], 'Aromatic'),
    **dict.fromkeys(['SER', 'THR', 'ASN', 'GLN', 'CYS'], 'Polar'),
    **dict.fromkeys(['LYS', 'ARG', 'HIS'], 'Positive'),
    **dict.fromkeys(['ASP', 'GLU'], 'Negative'),
    **dict.fromkeys(['DA', 'DT', 'DG', 'DC', 'A', 'U', 'G', 'C'], 'Nucleotide'),
}
# light tints of the palette, so residue boxes stay quiet behind the contact lines
CLASS_COLOR = {'Hydrophobic': '#EFEFEF', 'Aromatic': '#D8DDE7', 'Polar': '#CCECE7',
               'Positive': '#DBF1F7', 'Negative': '#FAD9D4', 'Nucleotide': '#FDEBE5', 'Other': '#FFFFFF'}
ELEMENT_COLOR = {'N': NAVY, 'O': RED, 'S': ORANGE, 'P': DARK, 'C': '#8C8C8C'}
CHAIN_COLORS = (BLUE, ORANGE)  # first chain, second chain: the same in every panel
CHAIN_CYCLE = (BLUE, RED, TEAL, NAVY, ORANGE, INK, DARK, GREY)


def _resnum(residue: str) -> int:
    m = re.search(r'-?\d+', residue)
    return int(m.group()) if m else 0


def _style(label: str):
    return TYPE_STYLE.get(label, (INK, '-', '?'))


def _short(res: str) -> str:
    name = res.split()[0]
    return f'{name.capitalize() if len(name) == 3 else name} {_resnum(res)}'


def _save(fig, output_file: str, **kwargs) -> str:
    """PNG at 300 dpi for reports, plus a vector PDF with editable text next to it."""
    fig.savefig(output_file, dpi=300, facecolor='white', **kwargs)
    fig.savefig(os.path.splitext(output_file)[0] + '.pdf', facecolor='white', **kwargs)
    plt.close(fig)
    return output_file


def _panel_x(fig, *axes) -> float:
    """Shared left edge (figure fraction) of a column of panels, tick labels included."""
    return min(ax.get_tightbbox(fig.canvas.get_renderer()).x0 for ax in axes) / fig.bbox.width


def _panel(fig, ax, letter: str, title: str, pad: float = 4, x: Optional[float] = None):
    """8 pt bold panel letter at the panel's true left edge, 7 pt heading beside it. Call after layout."""
    x = _panel_x(fig, ax) if x is None else x
    y = ax.get_position().y1 + pad / 72 / fig.get_figheight()
    fig.text(x, y, letter, fontsize=8, fontweight='bold', ha='left', va='bottom')
    fig.text(x + 12 / 72 / fig.get_figwidth(), y, title, fontsize=7, ha='left', va='bottom')


def _pair_table(df: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Unique (left residue, right residue, type) contacts between chains a and b."""
    fwd = df[(df['Chain1'] == a) & (df['Chain2'] == b)][['Residue1', 'Residue2', 'Interaction_Type', 'Distance_Å']]
    rev = df[(df['Chain1'] == b) & (df['Chain2'] == a)][['Residue2', 'Residue1', 'Interaction_Type', 'Distance_Å']]
    cols = ['left', 'right', 'type', 'dist']
    t = pd.concat([fwd.set_axis(cols, axis=1), rev.set_axis(cols, axis=1)])
    if a == b:  # intra-chain: same contact may be listed in both directions
        swap = t['left'].map(_resnum) > t['right'].map(_resnum)
        t.loc[swap, ['left', 'right']] = t.loc[swap, ['right', 'left']].values
    return t.groupby(['left', 'right', 'type'], as_index=False)['dist'].min()


def _contacts(df: pd.DataFrame, a: str, b: str, max_residues: int):
    """Contacts to draw for a chain pair: (table, note). vdW only when nothing else exists."""
    t = _pair_table(df, a, b)
    specific = t[t['type'] != 'Van der Waals']
    t = specific if len(specific) else t
    note = ''
    for side in ('left', 'right'):  # keep the busiest residues so the figure stays readable
        counts = t[side].value_counts()
        if len(counts) > max_residues:
            t = t[t[side].isin(counts.index[:max_residues])]
            note = f'the {max_residues} most connected residues per chain'
    return t, note


def _node(ax, x, y, res, chain, buried, ha='center'):
    deep = buried.get((chain, res), 0) >= BURIED_MIN
    ax.text(x, y, _short(res), ha=ha, va='center', fontsize=6, fontweight='bold', zorder=4,
            bbox=dict(boxstyle='round,pad=0.3', fc=CLASS_COLOR[RESIDUE_CLASS.get(res.split()[0], 'Other')],
                      ec=INK, lw=1.2 if deep else 0.4))


def _badge(ax, x, y, text, color):
    ax.text(x, y, text, fontsize=5, color=color, ha='center', va='center', zorder=3,
            bbox=dict(boxstyle='circle,pad=0.2', fc='white', ec=color, lw=0.5))


def _distance(ax, x, y, row, color, **kwargs):
    # distances only for the directional contacts, to keep the map readable
    if row.type in ('Hydrogen bond', 'Salt bridge / ionic', 'Disulfide') and not np.isnan(row.dist):
        ax.text(x, y, f'{row.dist:.1f} Å', fontsize=5, color=color, zorder=2,
                bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85), **kwargs)


def _draw_columns(ax, t: pd.DataFrame, chain_a: str, chain_b: str, buried: Dict):
    """Two-column residue contact map with a type badge on every line."""
    left = sorted(t['left'].unique(), key=_resnum)
    lpos = {r: i for i, r in enumerate(left)}
    # order the right side by where its partners sit, to cut line crossings
    right = sorted(t['right'].unique(), key=lambda r: (t[t['right'] == r]['left'].map(lpos).mean(), _resnum(r)))
    n = max(len(left), len(right))
    ly = {r: n - 1 - i * (n - 1) / max(len(left) - 1, 1) for i, r in enumerate(left)}
    ry = {r: n - 1 - i * (n - 1) / max(len(right) - 1, 1) for i, r in enumerate(right)}
    if len(left) == 1:
        ly[left[0]] = (n - 1) / 2
    if len(right) == 1:
        ry[right[0]] = (n - 1) / 2

    s = np.linspace(0, 1, 60)
    rank = list(TYPE_STYLE)
    for i, ((l, r), grp) in enumerate(t.groupby(['left', 'right'], sort=False)):
        grp = grp.sort_values('type', key=lambda c: c.map(lambda v: rank.index(v) if v in rank else 99))
        for k, row in enumerate(grp.itertuples()):
            bend = (k - (len(grp) - 1) / 2) * 0.5  # fan out several contact types between one pair
            y0, y1 = ly[l], ry[r]
            ys = (1 - s) ** 2 * y0 + 2 * (1 - s) * s * ((y0 + y1) / 2 + bend) + s ** 2 * y1
            color, ls, badge = _style(row.type)
            ax.plot(s, ys, color=color, ls=ls, lw=0.8, zorder=1)
            j = int(len(s) * (0.22 + 0.5 * ((i * 0.618 + k * 0.31) % 1)))  # stagger badges along the line
            _badge(ax, s[j], ys[j], badge, color)
            j2 = min(j + 4, len(s) - 1)
            _distance(ax, s[j2], ys[j2], row, color, ha='left', va='bottom')

    for chain, pos, x, ha in ((chain_a, ly, 0.0, 'right'), (chain_b, ry, 1.0, 'left')):
        for res, y in pos.items():
            _node(ax, x, y, res, chain, buried, ha=ha)
    ax.text(0.0, n + 0.2, f'Chain {chain_a}', ha='right', va='bottom', fontsize=7, fontweight='bold')
    ax.text(1.0, n + 0.2, f'Chain {chain_b}', ha='left', va='bottom', fontsize=7, fontweight='bold')
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-1, n + 1)
    ax.axis('off')
    return left + right


def _draw_radial(ax, t: pd.DataFrame, chain_a: str, chain_b: str, peptide: str, buried: Dict):
    """PandaMap-style map: the peptide as a backbone in the centre, partner residues around it."""
    pside, oside = ('left', 'right') if peptide == chain_a else ('right', 'left')
    other = chain_b if peptide == chain_a else chain_a
    # shortcut: only contacting peptide residues are drawn, pass the full sequence in if gaps matter
    pep = sorted(t[pside].unique(), key=_resnum)
    step = 2.0
    ppos = {r: ((i - (len(pep) - 1) / 2) * step, 0.5 if i % 2 else -0.5) for i, r in enumerate(pep)}
    px = {r: xy[0] for r, xy in ppos.items()}
    prot = sorted(t[oside].unique(), key=lambda r: (t[t[oside] == r][pside].map(px).mean(), _resnum(r)))
    rx = max(len(pep) * step / 2 + 3.0, len(prot) * 0.55, 6.0)
    ry = rx * 0.7
    opos = {}
    for group, sign in ((prot[0::2], 1), (prot[1::2], -1)):  # alternate above / below the peptide
        for i, r in enumerate(group):
            ang = np.radians(168 - 156 * (i + 0.5) / len(group))
            opos[r] = (rx * np.cos(ang), sign * ry * (0.35 + 0.65 * np.sin(ang)))

    xs, ys = zip(*[ppos[r] for r in pep])
    ax.plot(xs, ys, color=INK, lw=2.5, solid_capstyle='round', solid_joinstyle='round', zorder=1)
    nucleic = RESIDUE_CLASS.get(pep[0].split()[0]) == 'Nucleotide'
    ends = ("5'", "3'") if nucleic else ('N', 'C')
    ax.text(xs[0] - 1.7, ys[0], ends[0], fontsize=6, ha='center', va='center', fontweight='bold')
    ax.text(xs[-1] + 1.7, ys[-1], ends[1], fontsize=6, ha='center', va='center', fontweight='bold')

    rank = list(TYPE_STYLE)
    for (pr, orr), grp in t.groupby([pside, oside], sort=False):
        grp = grp.sort_values('type', key=lambda c: c.map(lambda v: rank.index(v) if v in rank else 99))
        (x0, y0), (x1, y1) = ppos[pr], opos[orr]
        for k, row in enumerate(grp.itertuples()):
            color, ls, badge = _style(row.type)
            ax.plot([x0, x1], [y0, y1], color=color, ls=ls, lw=0.8, zorder=2)
            f = 0.78 - 0.16 * k  # badges of one residue pair sit in a row along its line
            _badge(ax, x0 + f * (x1 - x0), y0 + f * (y1 - y0), badge, color)
            if k == 0:
                _distance(ax, x0 + 0.5 * (x1 - x0), y0 + 0.5 * (y1 - y0), row, color, ha='center', va='center')
    for r, (x, y) in ppos.items():
        _node(ax, x, y, r, peptide, buried)
    for r, (x, y) in opos.items():
        _node(ax, x, y, r, other, buried)

    ax.set_xlim(-rx - 1.8, rx + 1.8)
    ax.set_ylim(-ry - 1.2, ry + 1.8)
    ax.text(0, ry + 1.5, f"{'Oligonucleotide' if nucleic else 'Peptide'}: chain {peptide}  ·  partner: chain {other}",
            ha='center', va='center', fontsize=7, fontweight='bold')
    ax.axis('off')
    return pep + prot


def _draw_map(ax, t: pd.DataFrame, chain_a: str, chain_b: str, buried: Optional[Dict] = None,
              peptide: Optional[str] = None):
    """Contact map for a chain pair: radial when one chain is a peptide, two columns otherwise."""
    buried = buried or {}
    if peptide:
        residues = _draw_radial(ax, t, chain_a, chain_b, peptide, buried)
    else:
        residues = _draw_columns(ax, t, chain_a, chain_b, buried)
    used = [k for k in TYPE_STYLE if k in set(t['type'])] + sorted(set(t['type']) - set(TYPE_STYLE))
    classes = [c for c in CLASS_COLOR if c in {RESIDUE_CLASS.get(r.split()[0], 'Other') for r in residues}]
    leg = ax.legend(handles=[Line2D([0], [0], color=_style(k)[0], ls=_style(k)[1], lw=0.9,
                                    label=f'{_style(k)[2]}  {k}') for k in used],
                    title='Interaction (badge)', loc='upper center', bbox_to_anchor=(0.2, 0.02), ncol=2,
                    handlelength=2.2, columnspacing=1.0)
    ax.add_artist(leg)
    handles = [Patch(fc=CLASS_COLOR[c], ec=INK, lw=0.4, label=c) for c in classes]
    if buried:
        handles.append(Patch(fc='white', ec=INK, lw=1.2, label=f'Buries ≥ {BURIED_MIN:.0f} Å²'))
    ax.legend(handles=handles, title='Residue', loc='upper center', bbox_to_anchor=(0.86, 0.02), ncol=2,
              handlelength=1.2, columnspacing=1.0)


def _peptide_of(chain_a: str, chain_b: str, short_chains) -> Optional[str]:
    """The chain to draw as a peptide: only when exactly one of the pair is short."""
    hits = [c for c in (chain_a, chain_b) if c in (short_chains or ())]
    return hits[0] if len(hits) == 1 and chain_a != chain_b else None


def _map_height(t: pd.DataFrame, peptide: Optional[str]) -> float:
    """Inches the map needs so 6 pt residue labels never touch."""
    return 3.4 if peptide else max(2.4, 0.17 * max(t['left'].nunique(), t['right'].nunique()) + 0.8)


@nature
def draw_interface_map(df: pd.DataFrame, chain_a: str, chain_b: str, output_file: str,
                       title: Optional[str] = None, max_residues: int = 30,
                       buried: Optional[Dict] = None, short_chains=None) -> Optional[str]:
    """Draw the residue-residue contacts between two chains; returns the file or None if no contacts."""
    t, note = _contacts(df, chain_a, chain_b, max_residues)
    if t.empty:
        return None
    peptide = _peptide_of(chain_a, chain_b, short_chains)
    fig, ax = plt.subplots(figsize=(5.5 if peptide else 4.5, _map_height(t, peptide) + 0.9))
    _draw_map(ax, t, chain_a, chain_b, buried, peptide)
    ax.set_title(title or f'Interface map: chain {chain_a} - chain {chain_b} ({len(t)} contacts'
                 + (f', {note}' if note else '') + ')', fontsize=7, pad=6)
    return _save(fig, output_file, bbox_inches='tight', pad_inches=0.1)


def chain_pairs(df: pd.DataFrame) -> List[tuple]:
    """Chain pairs with specific (non-vdW) contacts, busiest first."""
    return list(df[df['Interaction_Type'] != 'Van der Waals'].apply(
        lambda r: tuple(sorted((r['Chain1'], r['Chain2']))), axis=1).value_counts().index)


@nature
def draw_interface_report(df: pd.DataFrame, chain_a: str, chain_b: str, output_file: str,
                          buried: Dict, name: str = '', max_residues: int = 25,
                          short_chains=None) -> Optional[str]:
    """One-page figure for a chain pair: a map, b contacts by type, c distances, d residue profile, e burial."""
    t, note = _contacts(df, chain_a, chain_b, max_residues)
    if t.empty:
        return None
    types = [k for k in TYPE_STYLE if k in set(t['type'])] + sorted(set(t['type']) - set(TYPE_STYLE))
    left = sorted(t['left'].unique(), key=_resnum)
    right = sorted(t['right'].unique(), key=_resnum)
    peptide = _peptide_of(chain_a, chain_b, short_chains)

    # layout in inches, so text keeps its point size whatever the number of residues
    map_h, low_h, head, gap, foot = _map_height(t, peptide), 1.9, 0.62, 1.15, 0.6
    height = head + map_h + gap + low_h + foot
    fig = plt.figure(figsize=(WIDTH, height))
    outer = fig.add_gridspec(2, 2, height_ratios=[map_h, low_h], width_ratios=[1.45, 1], wspace=0.5,
                             hspace=gap / ((map_h + low_h) / 2), left=0.15, right=0.98,
                             top=1 - head / height, bottom=foot / height)
    right_col = outer[0, 1].subgridspec(2, 1, hspace=0.75)

    # header: the numbers a reader wants first
    pair = df[((df['Chain1'] == chain_a) & (df['Chain2'] == chain_b)) |
              ((df['Chain1'] == chain_b) & (df['Chain2'] == chain_a))]
    bsa_a = sum(v for (c, _), v in buried.items() if c == chain_a)
    bsa_b = sum(v for (c, _), v in buried.items() if c == chain_b)
    facts = [f'{len(t)} contacts between {len(left)} + {len(right)} residues']
    hb = pair[pair['Interaction_Type'] == 'Hydrogen bond']
    if len(hb):
        c = hb.loc[hb['Distance_Å'].idxmin()]
        facts.append(f"{(t['type'] == 'Hydrogen bond').sum()} hydrogen bonds, closest {c['Residue1']} {c['Atom1']}-"
                     f"{c['Residue2']} {c['Atom2']} at {c['Distance_Å']:.2f} Å")
    if buried:
        facts.append(f'buried surface {bsa_a + bsa_b:,.0f} Å² ({chain_a}: {bsa_a:,.0f}, {chain_b}: {bsa_b:,.0f})')
    vdw = (pair['Interaction_Type'] == 'Van der Waals').sum()
    if vdw:
        facts.append(f'{vdw} van der Waals atom contacts (not drawn)')
    fig.text(0.02, 1 - 0.1 / height, f'{name} · interface chain {chain_a} - chain {chain_b}'.lstrip(' ·'),
             fontsize=7, fontweight='bold', va='top')
    fig.text(0.02, 1 - 0.23 / height, '  ·  '.join(facts), fontsize=5.5, va='top')

    ax_a = fig.add_subplot(outer[0, 0])
    _draw_map(ax_a, t, chain_a, chain_b, buried, peptide)

    ax_b = fig.add_subplot(right_col[0])
    counts = t['type'].value_counts().reindex(types)
    bars = ax_b.barh(types, counts.values, color=[_style(k)[0] for k in types], height=0.6)
    for bar, k in zip(bars, types):  # same line style as the map tells types of one colour apart
        bar.set(edgecolor=_style(k)[0], linewidth=0.5, alpha=1 if _style(k)[1] == '-' else 0.45)
    ax_b.bar_label(bars, padding=2, fontsize=5.5)
    ax_b.set_ylim(max(len(types), 5) - 0.5, -0.5)  # keeps bars slim when few types are present
    ax_b.set_xlabel('Residue-residue contacts (n)')

    ax_c = fig.add_subplot(right_col[1])
    dtypes = [k for k in types if t.loc[t['type'] == k, 'dist'].notna().any()]
    for i, k in enumerate(dtypes):
        d = np.sort(t.loc[t['type'] == k, 'dist'].dropna().values)
        ax_c.scatter(d, i + np.linspace(-0.22, 0.22, len(d)) if len(d) > 1 else [i], s=5, color=_style(k)[0],
                     linewidths=0, zorder=3)
        ax_c.plot([np.median(d)] * 2, [i - 0.34, i + 0.34], color=INK, lw=0.6)
    ax_c.set_yticks(range(len(dtypes)), dtypes)
    ax_c.set_ylim(max(len(dtypes), 5) - 0.5, -0.5)
    ax_c.set_xlabel('Distance (Å); bar, median')

    ax_d = fig.add_subplot(outer[1, 0])
    cols = [(chain_a, r) for r in left] + [(chain_b, r) for r in right]
    xs = {c: i + (1 if c[0] == chain_b and chain_a != chain_b else 0) for i, c in enumerate(cols)}  # gap between chains
    for yi, k in enumerate(types):
        sub = t[t['type'] == k]
        hit = [c for c in cols if c in {(chain_a, r) for r in sub['left']} | {(chain_b, r) for r in sub['right']}]
        ax_d.scatter([xs[c] for c in hit], [yi] * len(hit), s=9, color=_style(k)[0], linewidths=0, zorder=3)
    rows = list(types)
    if buried:
        b = np.array([buried.get(c, 0) for c in cols])
        ax_d.scatter([xs[c] for c in cols], [len(types)] * len(cols), s=1 + b * 0.25, facecolor='white',
                     edgecolor=INK, lw=0.5, zorder=3)
        rows.append('Buried area (circle size)')
    ax_d.set_yticks(range(len(rows)), rows)
    ax_d.set_xticks([xs[c] for c in cols], [_short(r) for _, r in cols], rotation=90, fontsize=5)
    ax_d.set_ylim(len(rows) - 0.4, -0.6)
    ax_d.set_xlim(-1, max(xs.values()) + 1)
    ax_d.grid(axis='y', color='#E5E5E5')
    ax_d.spines[['left', 'bottom']].set_visible(False)
    ax_d.tick_params(length=0)
    for side, (chain, members) in enumerate(((chain_a, left), (chain_b, right))):
        if members:
            mid = np.mean([xs[(chain, r)] for r in members])
            ax_d.text(mid, 1.02, f'Chain {chain}', ha='center', va='bottom', fontsize=6, fontweight='bold',
                      color=INK, transform=ax_d.get_xaxis_transform())

    ax_e = fig.add_subplot(outer[1, 1])
    if buried:
        top = sorted(buried.items(), key=lambda kv: kv[1], reverse=True)[:14][::-1]
        ax_e.barh([f'{c}:{_short(r)}' for (c, r), _ in top], [v for _, v in top],
                  color=[CHAIN_COLORS[0] if c == chain_a else CHAIN_COLORS[1] for (c, _), _ in top], height=0.65)
        ax_e.tick_params(axis='y', labelsize=5.5)
        ax_e.legend(handles=[Patch(fc=CHAIN_COLORS[0], label=f'Chain {chain_a}'),
                             Patch(fc=CHAIN_COLORS[1], label=f'Chain {chain_b}')], loc='lower right',
                    handlelength=1.2)
        ax_e.set_xlabel('Surface buried on binding (Å²)')

    fig.canvas.draw()
    x_left, x_right = _panel_x(fig, ax_a, ax_d), _panel_x(fig, ax_b, ax_c, ax_e)  # one edge per column
    _panel(fig, ax_a, 'a', 'Interaction map' + (f' ({note})' if note else ''), x=x_left)
    _panel(fig, ax_b, 'b', 'Contacts by type', x=x_right)
    _panel(fig, ax_c, 'c', 'Contact distances', x=x_right)
    _panel(fig, ax_d, 'd', 'Residue profile', pad=14, x=x_left)
    _panel(fig, ax_e, 'e', 'Most buried residues', x=x_right)
    fig.text(0.98, 0.08 / height, 'PandaPPI. Burial: solvent-accessible area of the isolated chain minus that in '
             'the pair (Shrake-Rupley, 1.4 Å probe).', ha='right', va='bottom', fontsize=5)
    _save(fig, output_file)
    print(f'Interface report saved to {output_file}')
    return output_file


@nature
def draw_contacts_3d(model, df: pd.DataFrame, chain_a: str, chain_b: str, output_file: str,
                     name: str = '', max_residues: int = 12, gap: float = 5.0) -> Optional[str]:
    """Atom-level view of an interface: real side-chain geometry projected to 2D, contacts between atoms.

    The view looks along the interface: chain_a on the left, chain_b on the right. The two chains are
    pulled apart by `gap` Å each so the contact lines are readable; printed distances are the real ones.
    """
    pair = df[(df['Interaction_Type'] != 'Van der Waals') &
              (((df['Chain1'] == chain_a) & (df['Chain2'] == chain_b)) |
               ((df['Chain1'] == chain_b) & (df['Chain2'] == chain_a)))]
    if pair.empty or chain_a == chain_b:
        return None
    flip = pair['Chain1'] == chain_b  # orient every row as chain_a -> chain_b
    c = pd.DataFrame({
        'ra': pair['Residue1'].where(~flip, pair['Residue2']), 'aa': pair['Atom1'].where(~flip, pair['Atom2']),
        'rb': pair['Residue2'].where(~flip, pair['Residue1']), 'ab': pair['Atom2'].where(~flip, pair['Atom1']),
        'type': pair['Interaction_Type'], 'dist': pair['Distance_Å']})
    rank = list(TYPE_STYLE)
    c['rank'] = c['type'].map(lambda v: rank.index(v) if v in rank else 99)
    # one line per residue pair: its most specific, then closest, contact
    c = c.sort_values(['rank', 'dist']).drop_duplicates(['ra', 'rb'])
    for col in ('ra', 'rb'):  # keep the busiest residues
        c = c[c[col].isin(c[col].value_counts().index[:max_residues])]

    # shortcut: residues are matched by name + number, insertion codes are not told apart
    lookup = {(ch.id, f"{r.resname} {r.id[1]}"): r for ch in model if ch.id in (chain_a, chain_b)
              for r in ch if r.id[0] == ' '}
    res_a = {r: lookup[(chain_a, r)] for r in c['ra'].unique() if (chain_a, r) in lookup}
    res_b = {r: lookup[(chain_b, r)] for r in c['rb'].unique() if (chain_b, r) in lookup}
    c = c[c['ra'].isin(res_a) & c['rb'].isin(res_b)]
    if c.empty:
        return None

    def heavy(res):
        return [a for a in res if a.element != 'H']

    xyz_a = np.array([a.coord for r in res_a.values() for a in heavy(r)])
    xyz_b = np.array([a.coord for r in res_b.values() for a in heavy(r)])
    mid = (xyz_a.mean(axis=0) + xyz_b.mean(axis=0)) / 2
    ex = xyz_b.mean(axis=0) - xyz_a.mean(axis=0)
    ex /= np.linalg.norm(ex)
    rest = np.vstack([xyz_a, xyz_b]) - mid
    rest = rest - np.outer(rest @ ex, ex)
    ey = np.linalg.svd(rest, full_matrices=False)[2][0]  # widest direction across the interface

    def project(coord, side):
        v = np.asarray(coord) - mid
        return np.array([v @ ex + (gap if side else -gap), v @ ey])

    fig, ax = plt.subplots(figsize=(WIDTH, 5.4))
    centers, extent = {}, [[], []]
    for side, residues in enumerate((res_a, res_b)):
        for rname, res in residues.items():
            atoms = heavy(res)
            pts = np.array([project(a.coord, side) for a in atoms])
            raw = np.array([a.coord for a in atoms])
            for i in range(len(atoms)):  # bonds guessed from distance, enough for standard residues
                for j in range(i + 1, len(atoms)):
                    if np.linalg.norm(raw[i] - raw[j]) < 1.9:
                        ax.plot(pts[[i, j], 0], pts[[i, j], 1], color=CHAIN_COLORS[side], lw=1.8, alpha=0.5,
                                solid_capstyle='round', zorder=1)
            centers[(side, rname)] = pts.mean(axis=0)
            extent[side].append(pts)

    def endpoint(side, rname, atom):
        res = (res_a, res_b)[side][rname]
        if atom and atom in res:
            return project(res[atom].coord, side), res[atom].element, atom
        ring = [a.coord for a in heavy(res) if a.get_id() not in ('N', 'CA', 'C', 'O')] or [a.coord for a in heavy(res)]
        return project(np.mean(ring, axis=0), side), None, ''  # ring / group contacts have no single atom

    for row in c.itertuples():
        color, ls, _ = _style(row.type)
        (p0, e0, n0), (p1, e1, n1) = endpoint(0, row.ra, row.aa), endpoint(1, row.rb, row.ab)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, ls=ls, lw=0.9, zorder=2)
        if not np.isnan(row.dist):
            m = p0 + 0.5 * (p1 - p0)
            ax.text(m[0], m[1], f'{row.dist:.1f} Å', fontsize=5, color=color, ha='center', va='center', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.9))
        for p, el, nm in ((p0, e0, n0), (p1, e1, n1)):
            ax.scatter(*p, s=16 if el else 9, color=ELEMENT_COLOR.get(el, '#8C8C8C'), edgecolor='white', lw=0.4,
                       zorder=4)
            if nm:
                ax.text(p[0], p[1] + 0.5, nm, fontsize=5, ha='center', va='bottom', zorder=5)

    # residue labels in a column at each outer edge, joined to their residue by a thin leader
    for side, (chain, residues) in enumerate(((chain_a, res_a), (chain_b, res_b))):
        allpts = np.vstack(extent[side])
        edge = allpts[:, 0].min() - 2.5 if side == 0 else allpts[:, 0].max() + 2.5
        names = sorted(residues, key=lambda r: centers[(side, r)][1])
        ys = [centers[(side, r)][1] for r in names]
        for i in range(1, len(ys)):  # push labels apart so they never overlap
            ys[i] = max(ys[i], ys[i - 1] + 1.7)
        shift = np.mean([centers[(side, r)][1] for r in names]) - np.mean(ys)
        for rname, y in zip(names, ys):
            cx, cy = centers[(side, rname)]
            ax.plot([edge, cx], [y + shift, cy], color=GREY, lw=0.3, zorder=0)
            ax.text(edge, y + shift, _short(rname), ha='right' if side == 0 else 'left', va='center', fontsize=6,
                    fontweight='bold', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.3', fc=CLASS_COLOR[RESIDUE_CLASS.get(rname.split()[0], 'Other')],
                              ec=CHAIN_COLORS[side], lw=0.8))
        ax.text(edge, max(ys) + shift + 2.4, f'Chain {chain}', ha='right' if side == 0 else 'left',
                fontsize=7, fontweight='bold')

    ax.set_aspect('equal')
    ax.margins(0.08)
    ax.axis('off')
    ax.set_title(f'{name} · atom-level contacts, chain {chain_a} - chain {chain_b}'.lstrip(' ·'), fontsize=7,
                 fontweight='bold', loc='left', pad=4)
    used = [k for k in TYPE_STYLE if k in set(c['type'])]
    handles = [Line2D([0], [0], color=_style(k)[0], ls=_style(k)[1], lw=0.9, label=k) for k in used]
    handles += [Line2D([0], [0], color=CHAIN_COLORS[0], lw=2, alpha=0.5, label=f'Chain {chain_a} residues'),
                Line2D([0], [0], color=CHAIN_COLORS[1], lw=2, alpha=0.5, label=f'Chain {chain_b} residues')]
    handles += [Line2D([0], [0], marker='o', ls='', mfc=ELEMENT_COLOR[e], mec='white', mew=0.4, ms=4, label=lab)
                for e, lab in (('N', 'Nitrogen'), ('O', 'Oxygen'), ('C', 'Carbon / ring centre'))]
    fig.legend(handles=handles, loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.035), handlelength=2.2)
    fig.text(0.5, 0.012, f'Atom positions projected along the interface. The chains are drawn {2 * gap:.0f} Å '
             'further apart than they are; printed distances are the measured ones.', ha='center', fontsize=5)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.16, top=0.95)
    _save(fig, output_file)
    print(f'3D contacts view saved to {output_file}')
    return output_file


def draw_all_interface_maps(df: pd.DataFrame, output_prefix: str, max_maps: int = 12,
                            short_chains=None) -> Dict[str, str]:
    """One map per chain pair that has contacts, busiest pairs first. Returns {'A-B': file}."""
    files = {}
    for a, b in chain_pairs(df)[:max_maps]:  # shortcut: capped at 12 maps, raise max_maps for huge assemblies
        out = draw_interface_map(df, a, b, f'{output_prefix}_map_{a}-{b}.png', short_chains=short_chains)
        if out:
            files[f'{a}-{b}'] = out
            print(f'Interface map saved to {out}')
    return files


@nature
def draw_dashboard(df: pd.DataFrame, output_file: str, title: str = 'PandaPPI interaction summary') -> str:
    """Four-panel overview: counts by type, chain-pair matrix, distance spread, hotspot residues."""
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.2))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.08, wspace=0.55, hspace=0.5)
    fig.text(0.02, 0.985, title, fontsize=7, fontweight='bold', va='top')

    counts = df['Interaction_Type'].value_counts().sort_values(ascending=False)
    ax = axes[0, 0]
    bars = ax.barh(counts.index, counts.values, color=[_style(k)[0] for k in counts.index], height=0.65)
    for bar, k in zip(bars, counts.index):
        bar.set(edgecolor=_style(k)[0], linewidth=0.5, alpha=1 if _style(k)[1] == '-' else 0.45)
    ax.bar_label(bars, padding=2, fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xscale('log')  # van der Waals contacts would otherwise flatten every other bar
    ax.set_xlabel('Contacts (n, log scale)')

    ax = axes[0, 1]
    chains = sorted(set(df['Chain1']) | set(df['Chain2']))
    m = pd.crosstab(df['Chain1'], df['Chain2']).reindex(index=chains, columns=chains, fill_value=0)
    m = m + m.T - np.diag(np.diag(m))  # each contact is stored once; make the matrix symmetric
    im = ax.imshow(m.values, cmap='Blues')
    ax.set_xticks(range(len(chains)), chains)
    ax.set_yticks(range(len(chains)), chains)
    ax.set_xlabel('Chain')
    ax.set_ylabel('Chain')
    ax.tick_params(length=0)
    for i in range(len(chains)):
        for j in range(len(chains)):
            v = m.values[i, j]
            if v:
                ax.text(j, i, int(v), ha='center', va='center', fontsize=5,
                        color='white' if v > m.values.max() / 2 else INK)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label('Contacts (n)', fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, width=0.5, length=2)
    cbar.outline.set_linewidth(0.5)

    ax = axes[1, 0]
    order = [k for k in counts.index if df.loc[df['Interaction_Type'] == k, 'Distance_Å'].notna().any()]
    line = dict(linewidth=0.5, color=INK)
    box = ax.boxplot([df.loc[df['Interaction_Type'] == k, 'Distance_Å'].dropna() for k in order], patch_artist=True,
                     vert=False, widths=0.6, boxprops=dict(line), whiskerprops=dict(line), capprops=dict(line),
                     medianprops=dict(linewidth=0.8, color=INK),
                     flierprops=dict(marker='o', markersize=1, markerfacecolor=GREY, markeredgewidth=0))
    ax.set_yticks(range(1, len(order) + 1), order)
    for patch, k in zip(box['boxes'], order):
        patch.set(facecolor=_style(k)[0], alpha=1 if _style(k)[1] == '-' else 0.45)
    ax.invert_yaxis()
    ax.set_xlabel('Distance (Å)')

    ax = axes[1, 1]
    spec = df[df['Interaction_Type'] != 'Van der Waals']
    spec = spec if len(spec) else df
    res = pd.concat([spec['Chain1'] + ':' + spec['Residue1'], spec['Chain2'] + ':' + spec['Residue2']])
    hot = res.value_counts().head(15)
    classes = [RESIDUE_CLASS.get(r.split(':')[1].split()[0], 'Other') for r in hot.index]
    ax.barh([f"{r.split(':')[0]}:{_short(r.split(':')[1])}" for r in hot.index], hot.values,
            color=[CLASS_COLOR[c] for c in classes], edgecolor=INK, lw=0.4, height=0.7)
    ax.invert_yaxis()
    ax.tick_params(axis='y', labelsize=5.5)
    ax.legend(handles=[Patch(fc=CLASS_COLOR[c], ec=INK, lw=0.4, label=c) for c in dict.fromkeys(classes)],
              loc='lower right', handlelength=1.2)
    ax.set_xlabel('Specific contacts (n)')

    fig.canvas.draw()
    edges = [_panel_x(fig, axes[0, i], axes[1, i]) for i in (0, 1)]  # one edge per column
    for ax, letter, heading in zip(axes.flat, 'abcd', ('Contacts by type', 'Contacts between chains',
                                                      'Contact distances', 'Hotspot residues')):
        _panel(fig, ax, letter, heading, x=edges[list(axes.flat).index(ax) % 2])
    _save(fig, output_file)
    print(f'Dashboard saved to {output_file}')
    return output_file
