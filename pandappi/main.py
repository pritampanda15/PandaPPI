# pandappi/main.py
"""
Command-line interface for PandaPPI.
"""

import argparse

from .core import PandaPPI

# CLI flag -> (interaction key, help text)
INTERACTION_FLAGS = {
    'hydrogen-bonds': ('hydrogen_bonds', 'hydrogen bonds'),
    'ionic': ('ionic_interactions', 'ionic interactions'),
    'hydrophobic': ('hydrophobic_interactions', 'hydrophobic interactions'),
    'pi-stacking': ('pi_stacking', 'pi-pi stacking interactions'),
    'pi-cation': ('pi_cation', 'pi-cation interactions'),
    'salt-bridges': ('salt_bridges', 'salt bridges'),
    'cation-pi': ('cation_pi', 'cation-pi interactions'),
    'ch-pi': ('ch_pi', 'CH-pi interactions'),
    'disulfide': ('disulfide_bridges', 'disulfide bridges'),
    'sulfur-aromatic': ('sulfur_aromatic', 'sulfur-aromatic interactions'),
    'water-mediated': ('water_mediated', 'water-mediated interactions'),
    'metal-coordination': ('metal_coordination', 'metal-coordinated bonds'),
    'halogen-bonds': ('halogen_bonds', 'halogen bonds'),
    'amide-aromatic': ('amide_aromatic', 'amide-aromatic interactions'),
    'van-der-waals': ('van_der_waals', 'van der Waals interactions'),
    'amide-amide': ('amide_amide', 'amide-amide hydrogen bonds'),
}
STANDARD = ['hydrogen_bonds', 'ionic_interactions', 'hydrophobic_interactions',
            'pi_stacking', 'pi_cation', 'salt_bridges']


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description='PandaPPI: Comprehensive Protein Interaction Mapper')
    parser.add_argument('pdb_file', help='Path to PDB file')

    parser.add_argument('--chains', nargs='+', help='Chains to analyze (e.g., A B)')
    parser.add_argument('--output', '-o', default='pandappi', help='Output file prefix (default: pandappi)')
    parser.add_argument('--all-outputs', action='store_true',
                        help='Generate everything: reports, interface maps, dashboard, network, 3D view, viewer scripts')

    out = parser.add_argument_group('Outputs')
    out.add_argument('--map', action='store_true', help='2D interface map (PNG) for every chain pair')
    out.add_argument('--interface-report', action='store_true',
                     help='Multi-panel report figure per chain pair: map, contact types, distances, '
                          'residue profile, buried surface area')
    out.add_argument('--contacts-3d', action='store_true', dest='contacts_3d',
                     help='Atom-level 3D contacts figure per chain pair (side chains, contacting atoms, distances)')
    out.add_argument('--dashboard', action='store_true', help='Summary dashboard figure (PNG)')
    out.add_argument('--html-report', action='store_true',
                     help='Full report: HTML page, text report, CSVs, dashboard and interface reports')
    out.add_argument('--report', action='store_true', help='Interaction table (CSV)')
    out.add_argument('--residue-summary', action='store_true', help='Per-residue interaction counts (CSV)')
    out.add_argument('--statistics', action='store_true', help='Print interaction statistics and save them as text')
    out.add_argument('--network', action='store_true', help='Residue interaction network (PNG + GML)')
    out.add_argument('--3d-plot', action='store_true', dest='plot_3d', help='Interactive 3D visualization (HTML)')
    out.add_argument('--export-vis', action='store_true',
                     help='Export visualization files for PyMOL, Chimera, VMD, and Molstar')

    parser.add_argument('--distance-cutoff', type=float, default=4.5,
                        help='Distance cutoff for interaction detection (default: 4.5Å)')
    parser.add_argument('--include-intrachain', action='store_true', help='Include intra-chain interactions')

    types = parser.add_argument_group('Interaction Types (default: all)')
    types.add_argument('--all-interactions', action='store_true', help='Map all interaction types (default)')
    types.add_argument('--standard-only', action='store_true',
                       help='Only standard interactions (H-bonds, ionic, hydrophobic, pi-stacking, pi-cation, salt bridges)')
    for flag, (_, text) in INTERACTION_FLAGS.items():
        types.add_argument(f'--{flag}', action='store_true', help=f'Map {text}')

    args = parser.parse_args()
    prefix = args.output

    analyzer = PandaPPI(args.pdb_file, args.chains)
    analyzer.map_interactions(distance_cutoff=args.distance_cutoff, include_intrachain=args.include_intrachain)

    selected = STANDARD if args.standard_only else [
        key for flag, (key, _) in INTERACTION_FLAGS.items() if getattr(args, flag.replace('-', '_'))]
    if selected and not args.all_interactions:
        print(f"Keeping selected interaction types: {', '.join(selected)}")
        analyzer.interactions = {k: v for k, v in analyzer.interactions.items() if k in selected}

    if not any(analyzer.interactions.values()):
        print("No interactions found - nothing to write.")
        return

    everything = args.all_outputs
    if args.html_report or everything:
        analyzer.generate_full_report(prefix)  # includes CSVs, interface reports and dashboard
    else:
        if args.report:
            analyzer.generate_report(f"{prefix}_report.csv")
        if args.residue_summary:
            from .reports import generator
            generator.export_to_csv(generator.create_residue_summary(analyzer.interactions),
                                    f"{prefix}_residue_summary.csv")
        if args.map:
            analyzer.create_interaction_maps(prefix)
        if args.interface_report:
            analyzer.create_interface_reports(prefix)
        if args.contacts_3d:
            analyzer.create_contact_views(prefix)
        if args.dashboard:
            analyzer.create_dashboard(f"{prefix}_dashboard.png")

    if args.network or everything:
        analyzer.create_interaction_network(f"{prefix}_network.png")
    if args.plot_3d or everything:
        analyzer.visualize_3d(f"{prefix}_3d.html")
    if args.export_vis or everything:
        vis_files = analyzer.export_visualization_scripts(prefix)
        print("\nExported visualization files:")
        for program, filename in vis_files.items():
            print(f" - {program}: {filename}")

    if args.statistics or everything:
        stats = analyzer.get_interaction_statistics()
        lines = ["Interaction Statistics:", f"Total interactions: {stats['total_interactions']}", "",
                 "Interactions by type:"]
        lines += [f"  - {k}: {v}" for k, v in stats['by_type'].items()]
        lines += ["", "Average distances:"]
        lines += [f"  - {k}: {v:.2f}Å" for k, v in stats['avg_distances'].items()]
        lines += ["", "Most interactive residues:"]
        residues = [f"  - {k}: {v} interactions" for k, v in stats['residue_frequencies'].items()]
        print("\n" + "\n".join(lines + residues[:10]))
        with open(f"{prefix}_stats.txt", 'w') as f:
            f.write("\n".join(lines + residues) + "\n")
        print(f"Statistics saved to {prefix}_stats.txt")

    print("Analysis complete!")


if __name__ == "__main__":
    main()
