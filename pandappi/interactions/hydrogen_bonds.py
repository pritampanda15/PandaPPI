"""
Module for detecting hydrogen bonds in protein structures.
"""

from typing import Dict, List, Tuple, Set, Optional
import math
import warnings

import numpy as np
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom 

from . import utils


def find_hydrogen_bonds(atoms_by_chain: Dict[str, List[Atom]], 
                       distance_cutoff: float = 3.5,
                       angle_cutoff: Optional[float] = None,
                       include_intrachain: bool = False) -> List[Dict]:
    """
    Find hydrogen bonds between atoms in different chains.
    
    Args:
        atoms_by_chain: Dictionary of atoms organized by chain
        distance_cutoff: Maximum distance for hydrogen bond detection (Å)
        angle_cutoff: Minimum angle for hydrogen bond detection (degrees)
        include_intrachain: Whether to include interactions within the same chain
        
    Returns:
        List of dictionaries containing hydrogen bond information
    """
    if angle_cutoff is not None:
        warnings.warn(
            "find_hydrogen_bonds: angle_cutoff is not applied -- detection is distance-based "
            "because hydrogens are usually absent from PDB files. The value is ignored.",
            UserWarning, stacklevel=2,
        )

    h_bonds = []
    
    # Define donors and acceptors
    # Donor atoms typically have covalently bound hydrogens
    donors = {
        'N': ['backbone', 'sidechain'],  # Backbone and sidechain nitrogens
        'O': ['sidechain'],            # Hydroxyl oxygens in sidechains
        'S': ['sidechain']             # Thiol groups
    }
    
    # Acceptor atoms typically have lone pairs
    acceptors = {
        'O': ['backbone', 'sidechain'],  # Carbonyl and hydroxyl oxygens
        'N': ['sidechain'],            # Some sidechain nitrogens
        'S': ['sidechain']             # Some sulfur atoms
    }
    
    # Process each UNORDERED pair of chains.
    # Iterating both (A, B) and (B, A) reported the same physical contact twice: `donors` and
    # `acceptors` both hold N, O and S, and only key membership is tested, so every heavy-atom
    # pair qualified in each direction.
    chain_ids = list(atoms_by_chain)
    chain_pairs = []
    for i, chain1 in enumerate(chain_ids):
        for chain2 in chain_ids[i:]:
            # Skip if same chain and not including intrachain interactions
            if chain1 == chain2 and not include_intrachain:
                continue

            chain_pairs.append((chain1, chain2))
    
    # Find hydrogen bonds
    for chain1, chain2 in chain_pairs:
        same_chain = chain1 == chain2
        for index1, atom1 in enumerate(atoms_by_chain[chain1]):
            # Check if atom1 is a potential donor
            if atom1.element in donors:
                # Within one chain, start past atom1: an atom is otherwise paired with itself
                # (distance 0.0, always under the cutoff) and every pair is seen twice.
                partners = (atoms_by_chain[chain2][index1 + 1:] if same_chain
                            else atoms_by_chain[chain2])
                # For each potential acceptor
                for atom2 in partners:
                    if atom2.element in acceptors:
                        # Calculate distance
                        dist = utils.calculate_distance(atom1, atom2)
                        
                        if dist <= distance_cutoff:
                            # For true hydrogen bonds, we'd need to find hydrogen atoms
                            # and calculate angles, but since H atoms are often missing
                            # in PDB files, we'll use a distance-based approximation
                            
                            atom1_res = atom1.get_parent()
                            atom2_res = atom2.get_parent()
                            
                            h_bonds.append({
                                'type': 'hydrogen_bond',
                                'donor_chain': chain1,
                                'donor_residue': f"{atom1_res.resname} {atom1_res.id[1]}",
                                'donor_atom': atom1.name,
                                'acceptor_chain': chain2,
                                'acceptor_residue': f"{atom2_res.resname} {atom2_res.id[1]}",
                                'acceptor_atom': atom2.name,
                                'distance': dist
                            })
    
    return h_bonds