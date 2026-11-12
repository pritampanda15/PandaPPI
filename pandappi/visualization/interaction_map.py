# pandappi/visualization/interaction_map.py
"""
2D figures for PandaPPI: PandaMap-style interface maps (one per chain pair)
and a summary dashboard. Both are drawn from the tidy interaction table
made by reports.generator.create_interaction_report().
"""

import re
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# label -> (colour, line style, badge); order is also drawing/legend priority
TYPE_STYLE = {
    'Hydrogen bond': ('#009e73', '-', 'H'),
    'Salt bridge / ionic': ('#e8472e', '-', 'I'),
    'Disulfide': ('#b8860b', '-', 'SS'),
    'Pi-pi stacking': ('#7b3fa0', '-', 'ππ'),
    'Cation-pi': ('#d6489a', '--', '+π'),
    'CH-pi': ('#8c564b', '--', 'Cπ'),
    'Hydrophobic': ('#8a8a8a', ':', 'h'),
    'Water-mediated': ('#0e7490', '-', 'W'),
    'Metal coordination': ('#3b1823', '-', 'M'),
    'Halogen bond': ('#00a087', '--', 'X'),
    'Sulfur-aromatic': ('#e08a00', '--', 'Sπ'),
    'Amide-aromatic': ('#4a78b0', '--', 'Aπ'),
    'Amide-amide': ('#6aa84f', '--', 'A'),
    'Van der Waals': ('#c7c7c7', ':', 'v'),
}
PEPTIDE_MAX = 30   # chains up to this many residues are drawn as a peptide (radial map)
BURIED_MIN = 40.0  # Å² a residue must bury to get the "buried" outline

RESIDUE_CLASS = {
    **dict.fromkeys(['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'GLY'], 'Hydrophobic'),
    **dict.fromkeys(['PHE', 'TYR', 'TRP'], 'Aromatic'),
    **dict.fromkeys(['SER', 'THR', 'ASN', 'GLN', 'CYS'], 'Polar'),
    **dict.fromkeys(['LYS', 'ARG', 'HIS'], 'Positive'),
    **dict.fromkeys(['ASP', 'GLU'], 'Negative'),
    **dict.fromkeys(['DA', 'DT', 'DG', 'DC', 'A', 'U', 'G', 'C'], 'Nucleotide'),
}
CLASS_COLOR = {'Hydrophobic': '#fde9a9', 'Aromatic': '#e3d0f5', 'Polar': '#cdeccd',
               'Positive': '#c9ddf7', 'Negative': '#f9cfcf', 'Nucleotide': '#ffd9b3', 'Other': '#e6e6e6'}


def _resnum(residue: str) -> int:
    m = re.search(r'-?\d+', residue)
    return int(m.group()) if m else 0


def _style(label: str):
    return TYPE_STYLE.get(label, ('#555555', '-', '?'))


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
            note = f'Showing the {max_residues} most connected residues per chain'
    return t, note


def _short(res: str) -> str:
    name = res.split()[0]
    return f'{name.capitalize() if len(name) == 3 else name} {_resnum(res)}'


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

    xl, xr = 0.0, 1.0
    s = np.linspace(0, 1, 60)
    rank = list(TYPE_STYLE)
    for i, ((l, r), grp) in enumerate(t.groupby(['left', 'right'], sort=False)):
        grp = grp.sort_values('type', key=lambda c: c.map(lambda v: rank.index(v) if v in rank else 99))
        for k, row in enumerate(grp.itertuples()):
            bend = (k - (len(grp) - 1) / 2) * 0.5  # fan out several contact types between one pair
            y0, y1 = ly[l], ry[r]
            ys = (1 - s) ** 2 * y0 + 2 * (1 - s) * s * ((y0 + y1) / 2 + bend) + s ** 2 * y1
            color, ls, badge = _style(row.type)
            ax.plot(xl + s * (xr - xl), ys, color=color, ls=ls, lw=1.5, alpha=0.9, zorder=1)
            j = int(len(s) * (0.22 + 0.5 * ((i * 0.618 + k * 0.31) % 1)))  # stagger badges along the line
            ax.text(s[j], ys[j], badge, fontsize=6.5, color=color, ha='center', va='center', zorder=3,
                    fontweight='bold', bbox=dict(boxstyle='circle,pad=0.25', fc='white', ec=color, lw=1))
            # distances only for the directional contacts, to keep the map readable
            if row.type in ('Hydrogen bond', 'Salt bridge / ionic', 'Disulfide') and not np.isnan(row.dist):
                j2 = min(j + 5, len(s) - 1)
                ax.text(s[j2], ys[j2], f'{row.dist:.1f} Å', fontsize=6, color=color, ha='left', va='bottom',
                        zorder=2, bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85))

    for chain, pos, x, ha in ((chain_a, ly, xl, 'right'), (chain_b, ry, xr, 'left')):
        for res, y in pos.items():
            _node(ax, x, y, res, chain, buried, ha=ha)

    ax.text(xl, n + 0.1, f'Chain {chain_a}', ha='right', va='bottom', fontsize=13, fontweight='bold')
    ax.text(xr, n + 0.1, f'Chain {chain_b}', ha='left', va='bottom', fontsize=13, fontweight='bold')
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-1, n + 1)
    ax.axis('off')
    return left + right



def _node(ax, x, y, res, chain, buried, ha='center', size=9):
    deep = buried.get((chain, res), 0) >= BURIED_MIN
    ax.text(x, y, _short(res), ha=ha, va='center', fontsize=size, fontweight='bold', zorder=4,
            bbox=dict(boxstyle='round,pad=0.35', fc=CLASS_COLOR[RESIDUE_CLASS.get(res.split()[0], 'Other')],
                      ec='#1f9bd1' if deep else '#555555', lw=2.4 if deep else 0.9))


def _draw_radial(ax, t: pd.DataFrame, chain_a: str, chain_b: str, peptide: str, buried: Dict):
    """PandaMap-style map: the peptide as a backbone in the centre, partner residues around it."""
    pside, oside = ('left', 'right') if peptide == chain_a else ('right', 'left')
    other = chain_b if peptide == chain_a else chain_a
    # shortcut: only contacting peptide residues are drawn, pass the full sequence in if gaps matter
    pep = sorted(t[pside].unique(), key=_resnum)
    step = 2.0
    ppos = {r: ((i - (len(pep) - 1) / 2) * step, 0.45 if i % 2 else -0.45) for i, r in enumerate(pep)}
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
    ax.plot(xs, ys, color='#33383d', lw=6, solid_capstyle='round', solid_joinstyle='round', zorder=1)
    nucleic = RESIDUE_CLASS.get(pep[0].split()[0]) == 'Nucleotide'
    ends = ("5'", "3'") if nucleic else ('N', 'C')
    ax.text(xs[0] - 1.5, ys[0], ends[0], fontsize=9, color='#33383d', ha='center', va='center', fontweight='bold')
    ax.text(xs[-1] + 1.5, ys[-1], ends[1], fontsize=9, color='#33383d', ha='center', va='center', fontweight='bold')

    rank = list(TYPE_STYLE)
    for (pr, orr), grp in t.groupby([pside, oside], sort=False):
        grp = grp.sort_values('type', key=lambda c: c.map(lambda v: rank.index(v) if v in rank else 99))
        (x0, y0), (x1, y1) = ppos[pr], opos[orr]
        for k, row in enumerate(grp.itertuples()):
            color, ls, badge = _style(row.type)
            ax.plot([x0, x1], [y0, y1], color=color, ls=ls, lw=1.5, alpha=0.9, zorder=2)
            f = 0.78 - 0.16 * k  # badges of one residue pair sit in a row along its line
            bx, by = x0 + f * (x1 - x0), y0 + f * (y1 - y0)
            ax.text(bx, by, badge, fontsize=6.5, color=color, ha='center', va='center', zorder=3,
                    fontweight='bold', bbox=dict(boxstyle='circle,pad=0.25', fc='white', ec=color, lw=1))
            if k == 0 and row.type in ('Hydrogen bond', 'Salt bridge / ionic', 'Disulfide') and not np.isnan(row.dist):
                ax.text(x0 + 0.5 * (x1 - x0), y0 + 0.5 * (y1 - y0), f'{row.dist:.1f} Å', fontsize=6, color=color,
                        ha='center', va='center', zorder=3,
                        bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85))
    for r, (x, y) in ppos.items():
        _node(ax, x, y, r, peptide, buried, size=9.5)
    for r, (x, y) in opos.items():
        _node(ax, x, y, r, other, buried)

    ax.set_xlim(-rx - 1.8, rx + 1.8)
    ax.set_ylim(-ry - 1.2, ry + 1.6)
    ax.text(0, ry + 1.3, f"{'Oligonucleotide' if nucleic else 'Peptide'}: chain {peptide}  ·  partner: chain {other}", ha='center', va='center',
            fontsize=12, fontweight='bold')
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
    leg1 = ax.legend(handles=[Line2D([0], [0], color=_style(k)[0], ls=_style(k)[1], lw=2,
                                     label=f'{_style(k)[2]}  {k}') for k in used],
                     title='Interaction (badge)', loc='upper center', bbox_to_anchor=(0.28, 0.02), ncol=2,
                     fontsize=8, title_fontsize=9, frameon=False)
    ax.add_artist(leg1)
    handles = [Patch(fc=CLASS_COLOR[c], ec='#555555', label=c) for c in classes]
    if buried:
        handles.append(Patch(fc='white', ec='#1f9bd1', lw=2.4, label=f'Buries ≥ {BURIED_MIN:.0f} Å²'))
    ax.legend(handles=handles, title='Residue', loc='upper center', bbox_to_anchor=(0.8, 0.02), ncol=2,
              fontsize=8, title_fontsize=9, frameon=False)




def _peptide_of(chain_a: str, chain_b: str, short_chains) -> Optional[str]:
    """The chain to draw as a peptide: only when exactly one of the pair is short."""
    hits = [c for c in (chain_a, chain_b) if c in (short_chains or ())]
    return hits[0] if len(hits) == 1 and chain_a != chain_b else None


def _map_rows(t: pd.DataFrame, peptide: Optional[str]) -> int:
    return 16 if peptide else max(t['left'].nunique(), t['right'].nunique())


def draw_interface_map(df: pd.DataFrame, chain_a: str, chain_b: str, output_file: str,
                       title: Optional[str] = None, max_residues: int = 30,
                       buried: Optional[Dict] = None, short_chains=None) -> Optional[str]:
    """Draw the residue-residue contacts between two chains; returns the file or None if no contacts."""
    t, note = _contacts(df, chain_a, chain_b, max_residues)
    if t.empty:
        return None
    peptide = _peptide_of(chain_a, chain_b, short_chains)
    n = _map_rows(t, peptide)
    fig, ax = plt.subplots(figsize=(12 if peptide else 10, max(5, 0.45 * n + 2.5)))
    _draw_map(ax, t, chain_a, chain_b, buried, peptide)
    ax.set_title(title or f'Interface map: chain {chain_a} - chain {chain_b}  ({len(t)} contacts)',
                 fontsize=13, pad=14)
    if note:
        fig.text(0.5, 0.005, note, ha='center', fontsize=8, style='italic', color='#555555')
    fig.savefig(output_file, dpi=200, bbox_inches='tight', pad_inches=0.4, facecolor='white')
    plt.close(fig)
    return output_file


def chain_pairs(df: pd.DataFrame) -> List[tuple]:
    """Chain pairs with specific (non-vdW) contacts, busiest first."""
    return list(df[df['Interaction_Type'] != 'Van der Waals'].apply(
        lambda r: tuple(sorted((r['Chain1'], r['Chain2']))), axis=1).value_counts().index)


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
    n = _map_rows(t, peptide)

    map_h = max(6.0, 0.42 * n + 2.5)
    fig = plt.figure(figsize=(17, map_h + 6.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[map_h / 2, map_h / 2, 4.6], width_ratios=[1.55, 1],
                          hspace=0.5, wspace=0.3, left=0.1, right=0.97, top=0.9, bottom=0.07)

    def panel(ax, letter, title, pad=10):
        ax.set_title(f'{letter}   {title}', loc='left', fontsize=13, fontweight='bold', pad=pad)
        ax.spines[['top', 'right']].set_visible(False)

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
    fig.text(0.04, 0.975, f'{name} · interface chain {chain_a} - chain {chain_b}'.lstrip(' ·'),
             fontsize=16, fontweight='bold', va='top')
    fig.text(0.04, 0.955, '   ·   '.join(facts), fontsize=10.5, va='top', color='#333333')

    ax = fig.add_subplot(gs[0:2, 0])
    _draw_map(ax, t, chain_a, chain_b, buried, peptide)
    panel(ax, 'a', 'Interaction map' + (f'  ({note.lower()})' if note else ''))

    ax = fig.add_subplot(gs[0, 1])
    counts = t['type'].value_counts().reindex(types)
    bars = ax.barh(types, counts.values, color=[_style(k)[0] for k in types], height=0.6)
    ax.bar_label(bars, padding=4, fontsize=9)
    ax.set_ylim(max(len(types), 5) - 0.5, -0.5)  # keeps bars slim when few types are present
    ax.set_xlabel('Residue-residue contacts (n)')
    panel(ax, 'b', 'Contacts by type')

    ax = fig.add_subplot(gs[1, 1])
    dtypes = [k for k in types if t.loc[t['type'] == k, 'dist'].notna().any()]
    for i, k in enumerate(dtypes):
        d = np.sort(t.loc[t['type'] == k, 'dist'].dropna().values)
        ax.scatter(d, i + np.linspace(-0.22, 0.22, len(d)) if len(d) > 1 else [i], s=18, color=_style(k)[0],
                   edgecolor='white', lw=0.4, zorder=3)
        ax.plot([np.median(d)] * 2, [i - 0.32, i + 0.32], color='#222222', lw=1)
    ax.set_yticks(range(len(dtypes)), dtypes)
    ax.invert_yaxis()
    ax.grid(axis='x', color='#e5e5e5')
    ax.set_xlabel('Distance (Å)   | = median')
    panel(ax, 'c', 'Contact distances')

    ax = fig.add_subplot(gs[2, 0])
    cols = [(chain_a, r) for r in left] + [(chain_b, r) for r in right]
    xs = {c: i + (1 if c[0] == chain_b and chain_a != chain_b else 0) for i, c in enumerate(cols)}  # gap between chains
    for yi, k in enumerate(types):
        sub = t[t['type'] == k]
        hit = {(chain_a, r) for r in sub['left']} | {(chain_b, r) for r in sub['right']}
        ax.scatter([xs[c] for c in cols if c in hit], [yi] * len([c for c in cols if c in hit]),
                   s=55, color=_style(k)[0], zorder=3)
    rows = list(types)
    if buried:
        b = np.array([buried.get(c, 0) for c in cols])
        ax.scatter([xs[c] for c in cols], [len(types)] * len(cols), s=4 + b * 1.6, facecolor='white',
                   edgecolor='#1f9bd1', lw=1.6, zorder=3)
        rows.append('Buried area (circle size)')
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xticks([xs[c] for c in cols], [_short(r) for _, r in cols], rotation=90, fontsize=8)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(-1, max(xs.values()) + 1)
    ax.grid(axis='y', color='#e5e5e5')
    ax.spines[['left', 'bottom']].set_visible(False)
    ax.tick_params(length=0)
    for chain, members in ((chain_a, left), (chain_b, right)):
        if members:
            mid = np.mean([xs[(chain, r)] for r in members])
            ax.text(mid, 1.02, f'Chain {chain}', ha='center', va='bottom', fontsize=10, fontweight='bold',
                    transform=ax.get_xaxis_transform())
    panel(ax, 'd', 'Residue profile', pad=28)

    ax = fig.add_subplot(gs[2, 1])
    if buried:
        top = sorted(buried.items(), key=lambda kv: kv[1], reverse=True)[:14][::-1]
        ax.barh([f'{c}:{_short(r)}' for (c, r), _ in top], [v for _, v in top],
                color=['#12233f' if c == chain_a else '#e08a00' for (c, _), _ in top], height=0.65)
        ax.legend(handles=[Patch(fc='#12233f', label=f'Chain {chain_a}'), Patch(fc='#e08a00', label=f'Chain {chain_b}')],
                  frameon=False, fontsize=9, loc='lower right')
        ax.set_xlabel('Surface buried on binding (Å²)')
    panel(ax, 'e', 'Most buried residues')

    fig.text(0.97, 0.012, 'PandaPPI · burial = solvent-accessible area of the isolated chain minus that in the pair '
             '(Shrake-Rupley, 1.4 Å probe)', ha='right', fontsize=8, color='#666666')
    fig.savefig(output_file, dpi=170, facecolor='white')
    plt.close(fig)
    print(f'Interface report saved to {output_file}')
    return output_file


ELEMENT_COLOR = {'N': '#2b4bd6', 'O': '#d62728', 'S': '#d4a800', 'P': '#e8650a', 'C': '#7f8790'}
CHAIN_COLORS = ('#3b6fb6', '#e08a00')  # left chain, right chain


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

    fig, ax = plt.subplots(figsize=(13, 11))
    centers, extent = {}, [[], []]
    for side, residues in enumerate((res_a, res_b)):
        for rname, res in residues.items():
            atoms = heavy(res)
            pts = np.array([project(a.coord, side) for a in atoms])
            raw = np.array([a.coord for a in atoms])
            for i in range(len(atoms)):  # bonds guessed from distance, enough for standard residues
                for j in range(i + 1, len(atoms)):
                    if np.linalg.norm(raw[i] - raw[j]) < 1.9:
                        ax.plot(pts[[i, j], 0], pts[[i, j], 1], color=CHAIN_COLORS[side], lw=4, alpha=0.38,
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
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, ls=ls, lw=1.8, zorder=2)
        if not np.isnan(row.dist):
            m = p0 + 0.5 * (p1 - p0)
            ax.text(m[0], m[1], f'{row.dist:.1f} Å', fontsize=7, color=color, ha='center', va='center', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.9))
        for (p, el, nm), dx in (((p0, e0, n0), -1), ((p1, e1, n1), 1)):
            ax.scatter(*p, s=110 if el else 60, color=ELEMENT_COLOR.get(el, '#7f8790'), edgecolor=color, lw=1.8,
                       zorder=4)
            if nm:
                ax.text(p[0], p[1] + 0.45, nm, fontsize=6.5, fontweight='bold', ha='center', va='bottom', zorder=5,
                        color='#222222')

    # residue labels in a column at each outer edge, joined to their residue by a thin leader
    for side, (chain, residues) in enumerate(((chain_a, res_a), (chain_b, res_b))):
        allpts = np.vstack(extent[side])
        edge = allpts[:, 0].min() - 2.5 if side == 0 else allpts[:, 0].max() + 2.5
        names = sorted(residues, key=lambda r: centers[(side, r)][1])
        ys = [centers[(side, r)][1] for r in names]
        for i in range(1, len(ys)):  # push labels apart so they never overlap
            ys[i] = max(ys[i], ys[i - 1] + 1.6)
        shift = (np.mean([centers[(side, r)][1] for r in names]) - np.mean(ys)) if names else 0
        for rname, y in zip(names, ys):
            cx, cy = centers[(side, rname)]
            ax.plot([edge, cx], [y + shift, cy], color='#b5b5b5', lw=0.7, zorder=0)
            ax.text(edge, y + shift, _short(rname), ha='right' if side == 0 else 'left', va='center', fontsize=9,
                    fontweight='bold', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.3', fc=CLASS_COLOR[RESIDUE_CLASS.get(rname.split()[0], 'Other')],
                              ec=CHAIN_COLORS[side], lw=1.6))
        ax.text(edge, max(ys) + shift + 2.2 if ys else 0, f'Chain {chain}', ha='right' if side == 0 else 'left',
                fontsize=13, fontweight='bold', color=CHAIN_COLORS[side])

    ax.set_aspect('equal')
    ax.margins(0.12)
    ax.axis('off')
    ax.set_title(f'3D contacts: {name} · chain {chain_a} - chain {chain_b}'.replace(':  ·', ':'), fontsize=15, pad=12)
    used = [k for k in TYPE_STYLE if k in set(c['type'])]
    handles = [Line2D([0], [0], color=_style(k)[0], ls=_style(k)[1], lw=2, label=k) for k in used]
    handles += [Line2D([0], [0], color=CHAIN_COLORS[0], lw=5, alpha=0.5, label=f'Chain {chain_a} residues'),
                Line2D([0], [0], color=CHAIN_COLORS[1], lw=5, alpha=0.5, label=f'Chain {chain_b} residues'),
                Line2D([0], [0], marker='o', ls='', mfc='#7f8790', mec='#333333', ms=9, label='Contacting atom / group')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.02))
    fig.text(0.5, 0.115, f'Real atom positions projected along the interface; the chains are pulled {2 * gap:.0f} Å '
             'apart for clarity, printed distances are the real ones.', ha='center', fontsize=8.5, color='#777777')
    fig.subplots_adjust(bottom=0.15, top=0.93)
    fig.savefig(output_file, dpi=170, facecolor='white', bbox_inches='tight', pad_inches=0.3)
    plt.close(fig)
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


def draw_dashboard(df: pd.DataFrame, output_file: str, title: str = 'PandaPPI interaction summary') -> str:
    """Four-panel overview: counts by type, chain-pair matrix, distance spread, hotspot residues."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle(title, fontsize=17, fontweight='bold')

    counts = df['Interaction_Type'].value_counts().sort_values()
    ax = axes[0, 0]
    bars = ax.barh(counts.index, counts.values, color=[_style(k)[0] for k in counts.index], edgecolor='#333333', lw=0.5)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_xscale('log')  # van der Waals contacts would otherwise flatten every other bar
    ax.set_title('Interactions by type (log scale)')
    ax.set_xlabel('Count')

    ax = axes[0, 1]
    chains = sorted(set(df['Chain1']) | set(df['Chain2']))
    m = pd.crosstab(df['Chain1'], df['Chain2']).reindex(index=chains, columns=chains, fill_value=0)
    m = m + m.T - np.diag(np.diag(m))  # each contact is stored once; make the matrix symmetric
    im = ax.imshow(m.values, cmap='YlGnBu')
    ax.set_xticks(range(len(chains)), chains)
    ax.set_yticks(range(len(chains)), chains)
    for i in range(len(chains)):
        for j in range(len(chains)):
            v = m.values[i, j]
            if v:
                ax.text(j, i, int(v), ha='center', va='center', fontsize=9,
                        color='white' if v > m.values.max() / 2 else 'black')
    ax.set_title('Contacts between chains')
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1, 0]
    order = [k for k in counts.index[::-1] if df.loc[df['Interaction_Type'] == k, 'Distance_Å'].notna().any()]
    box = ax.boxplot([df.loc[df['Interaction_Type'] == k, 'Distance_Å'].dropna() for k in order],
                     patch_artist=True, vert=False, flierprops=dict(markersize=2))
    ax.set_yticks(range(1, len(order) + 1), order)
    for patch, k in zip(box['boxes'], order):
        patch.set_facecolor(_style(k)[0])
        patch.set_alpha(0.75)
    ax.invert_yaxis()
    ax.set_title('Distance distribution')
    ax.set_xlabel('Distance (Å)')

    ax = axes[1, 1]
    spec = df[df['Interaction_Type'] != 'Van der Waals']
    spec = spec if len(spec) else df
    res = pd.concat([spec['Chain1'] + ':' + spec['Residue1'], spec['Chain2'] + ':' + spec['Residue2']])
    hot = res.value_counts().head(15).sort_values()
    colors = [CLASS_COLOR[RESIDUE_CLASS.get(r.split(':')[1].split()[0], 'Other')] for r in hot.index]
    ax.barh(hot.index, hot.values, color=colors, edgecolor='#333333', lw=0.6)
    ax.set_title('Hotspot residues (specific interactions)')
    ax.set_xlabel('Interactions')

    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_file, dpi=200, facecolor='white')
    plt.close(fig)
    print(f'Dashboard saved to {output_file}')
    return output_file
