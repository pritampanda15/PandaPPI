"""Unit checks for `find_hydrogen_bonds`: one contact is one row, and `angle_cutoff` is honest.

These use plain stub atoms rather than a PDB: the function only reads `.element`, `.name`,
`.coord` and `.get_parent()`, so a two-atom case pins the pairing logic exactly.
"""
import warnings

import numpy as np
import pytest

from pandappi.interactions import hydrogen_bonds


class _Residue:
    def __init__(self, resname, number):
        self.resname = resname
        self.id = (' ', number, ' ')


class _Atom:
    def __init__(self, name, element, coord, residue):
        self.name = name
        self.element = element
        self.coord = coord
        self._residue = residue

    def get_parent(self):
        return self._residue

    def __sub__(self, other):  # Bio.PDB.Atom gives the distance via subtraction
        return float(np.linalg.norm(np.subtract(self.coord, other.coord)))


def _two_chains(distance=3.0):
    """One N in chain A, one O in chain B, `distance` apart: a single physical contact."""
    asn, asp = _Residue('ASN', 33), _Residue('ASP', 99)
    return {
        'A': [_Atom('ND2', 'N', (0.0, 0.0, 0.0), asn)],
        'B': [_Atom('OD1', 'O', (distance, 0.0, 0.0), asp)],
    }


def test_one_contact_is_one_row():
    bonds = hydrogen_bonds.find_hydrogen_bonds(_two_chains())
    assert len(bonds) == 1, bonds  # both (A, B) and (B, A) used to be reported
    assert bonds[0]['donor_chain'] == 'A'
    assert bonds[0]['acceptor_chain'] == 'B'
    assert bonds[0]['distance'] == pytest.approx(3.0)


def test_distance_cutoff_is_applied():
    assert hydrogen_bonds.find_hydrogen_bonds(_two_chains(distance=3.4)) != []
    assert hydrogen_bonds.find_hydrogen_bonds(_two_chains(distance=3.6)) == []


def test_intrachain_pairs_each_atom_once_and_never_with_itself():
    asn, asp = _Residue('ASN', 33), _Residue('ASP', 99)
    atoms = {'A': [_Atom('ND2', 'N', (0.0, 0.0, 0.0), asn),
                   _Atom('OD1', 'O', (3.0, 0.0, 0.0), asp)]}
    bonds = hydrogen_bonds.find_hydrogen_bonds(atoms, include_intrachain=True)
    assert len(bonds) == 1, bonds              # was 4: two self-pairs plus the pair twice
    assert all(b['distance'] > 0.0 for b in bonds)


def test_angle_cutoff_warns_instead_of_being_silently_ignored():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        hydrogen_bonds.find_hydrogen_bonds(_two_chains(), angle_cutoff=150)
    assert any(issubclass(w.category, UserWarning) and 'angle_cutoff' in str(w.message)
               for w in caught)


def test_no_warning_by_default():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        hydrogen_bonds.find_hydrogen_bonds(_two_chains())
    assert [w for w in caught if issubclass(w.category, UserWarning)] == []
