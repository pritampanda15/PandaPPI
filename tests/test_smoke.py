"""End-to-end check on the smallest example: every report file is written and the numbers are sane."""
import os

from pandappi import PandaPPI

PDB = os.path.join(os.path.dirname(__file__), '..', 'test', 'protein-nucleicacid', '3l1p.pdb')


def test_full_report(tmp_path):
    analyzer = PandaPPI(PDB)
    analyzer.map_interactions()
    files = analyzer.generate_full_report(str(tmp_path / 'out'))
    for key in ('interactions_csv', 'residues_csv', 'text', 'dashboard', 'html', 'interface_A-N', 'contacts3d_A-N'):
        assert os.path.getsize(files[key]) > 0, key
    buried = analyzer.buried_surface('A', 'N')
    assert 500 < sum(buried.values()) < 5000  # a protein-DNA interface, not zero and not the whole surface
    assert analyzer._short_chains() >= {'M', 'N'}  # the DNA strands get the radial layout
