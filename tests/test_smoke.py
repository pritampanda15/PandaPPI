"""End-to-end check on the smallest example: every output is written and the numbers are sane."""
import os

import pytest

from pandappi import PandaPPI

PDB = os.path.join(os.path.dirname(__file__), '..', 'test', 'protein-nucleicacid', '3l1p.pdb')


@pytest.fixture(scope='module')
def analyzer():
    a = PandaPPI(PDB)
    a.map_interactions()
    return a


def test_full_report(analyzer, tmp_path):
    files = analyzer.generate_full_report(str(tmp_path / 'out'))
    for key in ('interactions_csv', 'residues_csv', 'text', 'dashboard', 'html', 'interface_A-N', 'contacts3d_A-N'):
        assert os.path.getsize(files[key]) > 0, key
    assert os.path.exists(files['interface_A-N'].replace('.png', '.pdf'))  # vector copy of every figure
    buried = analyzer.buried_surface('A', 'N')
    assert 500 < sum(buried.values()) < 5000  # a protein-DNA interface, not zero and not the whole surface
    assert analyzer._short_chains() >= {'M', 'N'}  # the DNA strands get the radial layout
    text = open(files['text'], encoding='utf-8').read()
    for heading in ('INTERACTIONS BY TYPE', 'INTERFACES', 'CLOSEST POLAR CONTACTS', 'HOTSPOT RESIDUES'):
        assert heading in text
    assert 'A:ARG 49   NH2   N:DG 10    O6    2.15' in text  # the closest hydrogen bond in 3L1P


def test_viewer_scripts(analyzer, tmp_path):
    files = analyzer.export_visualization_scripts(str(tmp_path / 'out'))
    assert set(files) == {'PyMOL', 'ChimeraX', 'Chimera', 'VMD', 'Molstar'}
    read = {k: open(v, encoding='utf-8').read() for k, v in files.items() if k != 'Molstar'}
    assert 'distance hydrogen_bond_001, (chain A and resi 41 and name O), (chain N and resi 10 and name OP1)' \
        in read['PyMOL']
    assert 'distance /A:41@O /N:10@OP1' in read['ChimeraX']
    assert 'distance :41.A@O :10.N@OP1' in read['Chimera']
    assert 'contact {chain A and resid "41" and name "O"} {chain N and resid "10" and name "OP1"}' in read['VMD']
    for text in read.values():
        assert 'out.pdb' in text  # each script loads the structure written next to it


def test_report_table_self_check():
    import runpy
    runpy.run_module('pandappi.reports.generator', run_name='__main__')
