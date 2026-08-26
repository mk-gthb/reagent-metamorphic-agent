"""Buried-position detection (relative SASA) + TM-score. Standard libs only — no Proto internals.
    pip install biotite tmtools numpy
"""
import io
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
from tmtools import tm_align

# Tien et al. 2013 theoretical max ASA (A^2) for relative-SASA normalization
MAX_ASA = {'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167, 'E': 223, 'Q': 225, 'G': 104,
           'H': 224, 'I': 197, 'L': 201, 'K': 236, 'M': 224, 'F': 240, 'P': 159, 'S': 155,
           'T': 172, 'W': 285, 'Y': 263, 'V': 174}
T2O = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q',
       'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
       'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def _load(text_or_path):
    """Accept raw PDB/CIF text or a path; auto-detect format; return amino-acid AtomArray (model 1)."""
    is_text = "\n" in text_or_path
    head = text_or_path[:400] if is_text else open(text_or_path).read(400)
    is_cif = ("_atom_site" in head) or head.lstrip().lower().startswith("data_")
    src = io.StringIO(text_or_path) if is_text else text_or_path
    if is_cif:
        arr = pdbx.get_structure(pdbx.CIFFile.read(src), model=1)
    else:
        arr = pdb.PDBFile.read(src).get_structure(model=1)
    return arr[struc.filter_amino_acids(arr)]


def buried_positions(pdb_in, rsasa_cutoff=0.15):
    """0-based residue indices (in structure order) with relative SASA < cutoff (buried core)."""
    arr = _load(pdb_in)
    atom_sasa = struc.sasa(arr, vdw_radii="Single")               # per-atom SASA (NaN for excluded)
    res_sasa = struc.apply_residue_wise(arr, atom_sasa, np.nansum)
    _, res_names = struc.get_residues(arr)
    rsasa = np.array([res_sasa[i] / MAX_ASA.get(T2O.get(n, 'A'), 129)
                      for i, n in enumerate(res_names)])
    return np.where(rsasa < rsasa_cutoff)[0].tolist()
    # CAVEAT: indices are structure-order; assumes 1 chain, no gaps -> aligns to seq string index.
    # If residue numbering has gaps, remap via arr.res_id before using as string offsets.


def _ca(pdb_in):
    arr = _load(pdb_in)
    ca = arr[arr.atom_name == "CA"]
    _, names = struc.get_residues(ca)
    return ca.coord, "".join(T2O.get(n, 'X') for n in names)


def tmscore(pdb_a, pdb_b):
    """TM-score between two structures (1.0 = identical fold)."""
    ca_a, sa = _ca(pdb_a)
    ca_b, sb = _ca(pdb_b)
    return tm_align(ca_a, ca_b, sa, sb).tm_norm_chain1
