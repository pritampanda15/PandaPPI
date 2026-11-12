# pandappi/visualization/network.py
"""
Residue interaction network for PandaPPI.
"""

from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D

from ..reports.generator import create_interaction_report
from .interaction_map import CHAIN_CYCLE, INK, TYPE_STYLE, WIDTH, _style, nature


@nature
def create_interaction_network(structure, interactions: Dict[str, List[Dict]]) -> Tuple[nx.Graph, plt.Figure]:
    """
    Build a residue-level graph (all interaction types) and draw it.

    Returns:
        Tuple of network graph and matplotlib figure
    """
    G = nx.Graph()
    df = create_interaction_report(interactions)
    if df is not None:
        for r in df.itertuples():
            n1, n2 = f"{r.Chain1}:{r.Residue1}", f"{r.Chain2}:{r.Residue2}"
            G.add_node(n1, chain=r.Chain1, residue=r.Residue1)
            G.add_node(n2, chain=r.Chain2, residue=r.Residue2)
            if G.has_edge(n1, n2):
                e = G[n1][n2]
                e['count'] += 1
                if r.Interaction_Type not in e['types']:
                    e['types'].append(r.Interaction_Type)
            else:
                G.add_edge(n1, n2, types=[r.Interaction_Type], count=1, distance=float(r[8]))

    fig, ax = plt.subplots(figsize=(WIDTH, 5.6))
    if G.number_of_nodes():
        rank = list(TYPE_STYLE)
        # an edge takes the colour of its most specific interaction type
        top = {(u, v): min(d['types'], key=lambda t: rank.index(t) if t in rank else 99)
               for u, v, d in G.edges(data=True)}
        chains = sorted({d['chain'] for _, d in G.nodes(data=True)})
        chain_color = {c: CHAIN_CYCLE[i % len(CHAIN_CYCLE)] for i, c in enumerate(chains)}
        pos = nx.spring_layout(G, k=1.5 / max(G.number_of_nodes(), 1) ** 0.5, iterations=80, seed=7)
        deg = dict(G.degree())

        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=[_style(top[e])[0] for e in G.edges()],
                               width=[0.3 if top[e] == 'Van der Waals' else 0.8 for e in G.edges()],
                               style=[_style(top[e])[1] for e in G.edges()])
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[chain_color[G.nodes[n]['chain']] for n in G],
                               node_size=[6 + 2.5 * deg[n] for n in G], edgecolors='white', linewidths=0.3)
        hubs = sorted(deg, key=deg.get, reverse=True)[:15]  # label only the hubs; the rest is in the GML
        nx.draw_networkx_labels(G, pos, ax=ax, labels={n: n for n in hubs}, font_size=5, font_color=INK)

        used = [t for t in rank if t in set(top.values())]
        leg = ax.legend(handles=[Line2D([0], [0], marker='o', ls='', color=chain_color[c], ms=3,
                                        label=f'Chain {c}') for c in chains], title='Chain', loc='upper left')
        ax.add_artist(leg)
        ax.legend(handles=[Line2D([0], [0], color=_style(t)[0], ls=_style(t)[1], lw=0.9, label=t) for t in used],
                  title='Interaction', loc='upper right')
    ax.set_title(f"Residue interaction network ({G.number_of_nodes()} residues, {G.number_of_edges()} contacts)",
                 fontsize=7, fontweight='bold', loc='left')
    ax.axis('off')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)
    return G, fig
