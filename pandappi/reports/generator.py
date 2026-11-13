# pandappi/reports/generator.py
"""
Reports for PandaPPI: one tidy table of every interaction, plus residue
summaries, a text report and a self-contained HTML report.
"""

import base64
import html
import os
import re
from typing import Dict, List, Optional

import pandas as pd

# Display name for each interaction key. Keys sharing a label are the same
# physical contact found by two detectors and are merged in diagrams.
LABELS = {
    'hydrogen_bonds': 'Hydrogen bond',
    'ionic_interactions': 'Salt bridge / ionic',
    'salt_bridges': 'Salt bridge / ionic',
    'hydrophobic_interactions': 'Hydrophobic',
    'pi_stacking': 'Pi-pi stacking',
    'pi_cation': 'Cation-pi',
    'cation_pi': 'Cation-pi',
    'ch_pi': 'CH-pi',
    'disulfide_bridges': 'Disulfide',
    'sulfur_aromatic': 'Sulfur-aromatic',
    'water_mediated': 'Water-mediated',
    'metal_coordination': 'Metal coordination',
    'halogen_bonds': 'Halogen bond',
    'amide_aromatic': 'Amide-aromatic',
    'van_der_waals': 'Van der Waals',
    'amide_amide': 'Amide-amide',
}

COLUMNS = ['Interaction_Type', 'Chain1', 'Residue1', 'Atom1',
           'Chain2', 'Residue2', 'Atom2', 'Distance_Å']

_CHAIN_KEY = re.compile(r'^(\w+_)?chain(\d?)$')


def _partners(interaction: Dict) -> List[tuple]:
    """(chain, residue, atom) for each side of an interaction dict.

    Every detector names its keys differently (donor_chain, pi_chain,
    chain1 ...), so the sides are found from the '*chain*' keys themselves.
    """
    sides = []
    for key, chain in interaction.items():
        m = _CHAIN_KEY.match(key)
        if not m or key.startswith('water'):
            continue
        prefix, num = m.group(1) or '', m.group(2)
        residue = interaction.get(f'{prefix}residue{num}')
        if residue is None:
            continue
        sides.append((chain, str(residue), str(interaction.get(f'{prefix}atom{num}', ''))))
    return sides


def create_interaction_report(interactions: Dict[str, List[Dict]]) -> Optional[pd.DataFrame]:
    """One row per detected interaction, for every interaction type."""
    rows, seen = [], set()
    for itype, items in interactions.items():
        for it in items:
            sides = _partners(it)
            if len(sides) < 2:  # e.g. a single residue coordinating a metal
                continue
            (c1, r1, a1), (c2, r2, a2) = sides[:2]
            label, dist = LABELS.get(itype, itype), round(float(it.get('distance', float('nan'))), 2)
            # two detectors (e.g. pi_cation / cation_pi) report the same contact: keep one
            key = (label, frozenset([(c1, r1, a1), (c2, r2, a2)]), dist)
            if key in seen:
                continue
            seen.add(key)
            rows.append((label, c1, r1, a1, c2, r2, a2, dist))
    if not rows:
        print("No interactions found to report.")
        return None
    return pd.DataFrame(rows, columns=COLUMNS)


def create_residue_summary(interactions: Dict[str, List[Dict]]) -> Optional[pd.DataFrame]:
    """Interaction counts per residue, one column per interaction type."""
    df = create_interaction_report(interactions)
    if df is None:
        return None
    both = pd.concat([
        df[['Chain1', 'Residue1', 'Interaction_Type']].set_axis(['Chain', 'Residue', 'Type'], axis=1),
        df[['Chain2', 'Residue2', 'Interaction_Type']].set_axis(['Chain', 'Residue', 'Type'], axis=1),
    ])
    summary = pd.crosstab([both['Chain'], both['Residue']], both['Type'])
    summary['Total'] = summary.sum(axis=1)
    return summary.sort_values('Total', ascending=False).reset_index()


def export_to_csv(df: pd.DataFrame, output_file: str) -> None:
    df.to_csv(output_file, index=False)
    print(f"Report exported to {output_file}")


def interface_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Interaction counts per chain pair (unordered), most contacts first."""
    pair = df.apply(lambda r: '-'.join(sorted((r['Chain1'], r['Chain2']))), axis=1)
    return pd.crosstab(pair, df['Interaction_Type']).assign(
        Total=lambda t: t.sum(axis=1)).sort_values('Total', ascending=False)


def _table(headers: list, rows: list, left: int = 1) -> str:
    """Aligned text table; the first `left` columns are left-aligned, the rest right-aligned."""
    rows = [[str(v) for v in r] for r in rows]
    widths = [max(len(v) for v in col) for col in zip(headers, *rows)]

    def line(r):
        return '  '.join(v.ljust(w) if i < left else v.rjust(w) for i, (v, w) in enumerate(zip(r, widths))).rstrip()
    return '\n'.join(['  ' + line(headers), '  ' + '  '.join('─' * w for w in widths)] + ['  ' + line(r) for r in rows])


def write_text_report(df: pd.DataFrame, summary: pd.DataFrame, pdb_file: str, output_file: str,
                      buried: Optional[Dict[str, float]] = None) -> None:
    """Plain-text report: summary, types, interfaces, closest contacts, hotspot residues.

    buried: optional {'A-B': Å² buried at that interface}.
    """
    from .. import __version__
    vdw = df['Interaction_Type'] == 'Van der Waals'
    spec = df[~vdw]
    pairs = interface_pairs(df)
    chains = sorted(set(df['Chain1']) | set(df['Chain2']))

    def section(title):
        return f"\n\n{title}\n{'━' * 72}\n"

    out = [f"PandaPPI {__version__}  ·  interaction report", '━' * 72,
           f"  Structure    {os.path.basename(pdb_file)}",
           f"  Chains       {', '.join(chains)}",
           f"  Contacts     {len(spec):,} specific  +  {int(vdw.sum()):,} van der Waals",
           f"  Residues     {len(summary):,} at an interface",
           f"  Interfaces   {len(pairs)} chain pairs"]
    if buried:
        out.append(f"  Buried area  {sum(buried.values()):,.0f} Å² in total")

    stats = df.groupby('Interaction_Type')['Distance_Å'].agg(['size', 'mean', 'min', 'max']).sort_values(
        'size', ascending=False)
    top = max(int(stats.loc[stats.index != 'Van der Waals', 'size'].max()) if len(spec) else 1, 1)
    rows = []
    for name, r in stats.iterrows():
        is_vdw = name == 'Van der Waals'
        num = ['–' if pd.isna(v) else f'{v:.2f}' for v in (r['mean'], r['min'], r['max'])]
        rows.append([name, f"{int(r['size']):,}", '–' if is_vdw else f"{100 * r['size'] / len(spec):.1f}%", *num,
                     ('' if is_vdw else '█' * max(1, round(20 * r['size'] / top))).ljust(20)])
    out.append(section('INTERACTIONS BY TYPE') + _table(
        ['Type', 'Count', 'Share', 'Mean Å', 'Min Å', 'Max Å', ''], rows)
        + "\n\n  Share is of specific contacts. Each contact is counted once.")

    def count(pair, col):
        return int(pairs.loc[pair, col]) if col in pairs.columns else 0
    rows = [[p, count(p, 'Total') - count(p, 'Van der Waals'), count(p, 'Hydrogen bond'),
             count(p, 'Salt bridge / ionic'), f"{count(p, 'Van der Waals'):,}",
             f"{buried[p]:,.0f}" if buried and p in buried else '–'] for p in pairs.index]
    out.append(section('INTERFACES') + _table(
        ['Chains', 'Specific', 'H-bonds', 'Salt bridges', 'van der Waals', 'Buried Å²'], rows))

    close = spec[spec['Interaction_Type'].isin(['Hydrogen bond', 'Salt bridge / ionic', 'Disulfide'])].nsmallest(
        12, 'Distance_Å')
    if len(close):
        rows = [[r.Interaction_Type, f"{r.Chain1}:{r.Residue1}", r.Atom1, f"{r.Chain2}:{r.Residue2}", r.Atom2,
                 f"{r[7]:.2f}"] for r in close.itertuples(index=False)]
        out.append(section('CLOSEST POLAR CONTACTS') + _table(
            ['Type', 'Residue 1', 'Atom', 'Residue 2', 'Atom', 'Å'], rows, left=5))

    hot = summary.copy()
    for col in ('Van der Waals', 'Hydrogen bond'):
        if col not in hot:
            hot[col] = 0
    hot['Specific'] = hot['Total'] - hot['Van der Waals']
    hot = hot.sort_values(['Specific', 'Total'], ascending=False).head(15)
    rows = [[f"{r['Chain']}:{r['Residue']}", r['Specific'], r['Hydrogen bond'], r['Van der Waals']]
            for _, r in hot.iterrows()]
    out.append(section('HOTSPOT RESIDUES') + _table(['Residue', 'Specific', 'H-bonds', 'van der Waals'], rows))

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out) + '\n')
    print(f"Text report saved to {output_file}")


_CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f5f6f8;color:#1c2330}
header{background:#12233f;color:#fff;padding:28px 40px}header h1{margin:0;font-size:26px}
header p{margin:6px 0 0;opacity:.8}main{max-width:1200px;margin:0 auto;padding:24px 40px}
.cards{display:flex;gap:16px;flex-wrap:wrap}.card{background:#fff;border-radius:10px;padding:16px 22px;
flex:1;min-width:150px;box-shadow:0 1px 3px rgba(0,0,0,.08)}.card b{display:block;font-size:28px}
section{background:#fff;border-radius:10px;padding:20px 24px;margin-top:20px;box-shadow:0 1px 3px rgba(0,0,0,.08)}
h2{margin-top:0;font-size:18px}img{max-width:100%;height:auto}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:5px 10px;text-align:left;
border-bottom:1px solid #e6e8ec}th{background:#f0f2f5;position:sticky;top:0}
.scroll{max-height:480px;overflow:auto}
"""


def write_html_report(df: pd.DataFrame, summary: pd.DataFrame, pdb_file: str,
                      images: Dict[str, str], output_file: str) -> None:
    """Self-contained HTML report; images (title -> PNG path) are embedded."""
    def img(path):
        with open(path, 'rb') as fh:
            return base64.b64encode(fh.read()).decode()

    def table(t, **kw):
        return f'<div class="scroll">{t.to_html(border=0, **kw)}</div>'

    name = html.escape(os.path.basename(pdb_file))
    cards = [('Interactions', len(df)), ('Interaction types', df['Interaction_Type'].nunique()),
             ('Residues involved', len(summary)), ('Chain pairs', len(interface_pairs(df)))]
    parts = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>PandaPPI report - {name}</title>',
             f'<style>{_CSS}</style></head><body><header><h1>PandaPPI interaction report</h1>',
             f'<p>{name}</p></header><main><div class="cards">']
    parts += [f'<div class="card"><b>{v}</b>{k}</div>' for k, v in cards]
    parts.append('</div>')
    for title, path in images.items():
        parts.append(f'<section><h2>{html.escape(title)}</h2>'
                     f'<img alt="{html.escape(title)}" src="data:image/png;base64,{img(path)}"></section>')
    parts.append('<section><h2>Interactions by chain pair</h2>' + table(interface_pairs(df)) + '</section>')
    parts.append('<section><h2>Hotspot residues</h2>' + table(summary.head(50), index=False) + '</section>')
    # shortcut: van der Waals rows are left out of the HTML table (thousands of rows), the CSV has them all
    detail = df[df['Interaction_Type'] != 'Van der Waals']
    parts.append('<section><h2>All specific interactions (van der Waals in CSV)</h2>'
                 + table(detail, index=False) + '</section></main></body></html>')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(''.join(parts))
    print(f"HTML report saved to {output_file}")


if __name__ == '__main__':  # self-check
    demo = {
        'hydrogen_bonds': [{'donor_chain': 'A', 'donor_residue': 'SER 1', 'donor_atom': 'OG',
                            'acceptor_chain': 'B', 'acceptor_residue': 'ASP 9', 'acceptor_atom': 'OD1',
                            'distance': 2.9}],
        'cation_pi': [{'cation_chain': 'A', 'cation_residue': 'LYS 3', 'cation_atom': 'NZ',
                       'pi_chain': 'B', 'pi_residue': 'TYR 7', 'distance': 4.1}],
        'water_mediated': [{'chain1': 'A', 'residue1': 'SER 1', 'chain2': 'B', 'residue2': 'ASP 9',
                            'water_chain': 'A', 'water_residue': 301, 'distance': 3.0}],
        'metal_coordination': [{'chain': 'A', 'residue': 'HIS 5', 'distance': 2.1}],
    }
    demo['pi_cation'] = [{'aromatic_chain': 'B', 'aromatic_residue': 'TYR 7', 'cationic_chain': 'A',
                          'cationic_residue': 'LYS 3', 'cationic_atom': 'NZ', 'distance': 4.1}]
    d = create_interaction_report(demo)
    assert len(d) == 3 and list(d.loc[1, ['Chain1', 'Residue1', 'Atom1', 'Chain2', 'Residue2']]) == \
        ['A', 'LYS 3', 'NZ', 'B', 'TYR 7'], d
    s = create_residue_summary(demo)
    assert s.loc[0, 'Total'] == 2 and interface_pairs(d).loc['A-B', 'Total'] == 3
    print('ok')
