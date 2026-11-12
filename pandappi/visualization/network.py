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
from .interaction_map import TYPE_STYLE, _style


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

    fig, ax = plt.subplots(figsize=(14, 11))
    if G.number_of_nodes():
        rank = list(TYPE_STYLE)
        # an edge takes the colour of its most specific interaction type
        top = {(u, v): min(d['types'], key=lambda t: rank.index(t) if t in rank else 99)
               for u, v, d in G.edges(data=True)}
        chains = sorted({d['chain'] for _, d in G.nodes(data=True)})
        cmap = plt.get_cmap('tab10')
        chain_color = {c: cmap(i % 10) for i, c in enumerate(chains)}
        pos = nx.spring_layout(G, k=1.5 / max(G.number_of_nodes(), 1) ** 0.5, iterations=80, seed=7)
        deg = dict(G.degree())

        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=[_style(top[e])[0] for e in G.edges()],
                               width=[0.6 if top[e] == 'Van der Waals' else 1.6 for e in G.edges()], alpha=0.65)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[chain_color[G.nodes[n]['chain']] for n in G],
                               node_size=[40 + 14 * deg[n] for n in G], edgecolors='white', linewidths=0.6)
        hubs = sorted(deg, key=deg.get, reverse=True)[:30]  # label only the hubs; the rest is in the GML
        nx.draw_networkx_labels(G, pos, ax=ax, labels={n: n for n in hubs}, font_size=7, font_weight='bold')

        used = [t for t in rank if t in set(top.values())]
        leg = ax.legend(handles=[Line2D([0], [0], marker='o', ls='', color=chain_color[c], label=f'Chain {c}')
                                 for c in chains], title='Chains', loc='upper left', frameon=False)
        ax.add_artist(leg)
        ax.legend(handles=[Line2D([0], [0], color=_style(t)[0], lw=2, label=t) for t in used],
                  title='Interaction', loc='upper right', frameon=False, fontsize=8)
    ax.set_title(f"Residue interaction network ({G.number_of_nodes()} residues, {G.number_of_edges()} contacts)",
                 fontsize=14)
    ax.axis('off')
    return G, fig
