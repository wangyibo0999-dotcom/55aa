"""Shared atom parsing and cysteine geometry; not a pipeline entry point."""
from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
import math
import numpy as np

@dataclass(frozen=True)
class Atom:
    name: str
    resname: str
    chain: str
    seq: int
    xyz: np.ndarray
    plddt: float

def parse_cif(path: Path) -> dict[tuple[str, int], list[Atom]]:
    residues: dict[tuple[str, int], list[Atom]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM "):
                continue
            fields = line.split()
            atom = Atom(
                name=fields[2],
                resname=fields[4],
                chain=fields[5],
                seq=int(fields[7]),
                plddt=float(fields[13]),
                xyz=np.array([float(fields[14]), float(fields[15]), float(fields[16])]),
            )
            residues[(atom.chain, atom.seq)].append(atom)
    return dict(residues)

def atom_map(atoms: list[Atom]) -> dict[str, Atom]:
    return {atom.name: atom for atom in atoms}

def norm(vector: np.ndarray) -> np.ndarray:
    length = np.linalg.norm(vector)
    if length == 0:
        raise ValueError("Zero-length vector")
    return vector / length

def place_sg(atoms: list[Atom], chi1_deg: float) -> np.ndarray:
    """Place SG using N-CA-CB-SG internal coordinates."""
    amap = atom_map(atoms)
    a = amap["N"].xyz
    b = amap["CA"].xyz
    c = amap["CB"].xyz
    bond_length = 1.81
    bond_angle = math.radians(114.0)
    phi = math.radians(chi1_deg)
    e1 = norm(b - c)
    plane_normal = norm(np.cross(b - a, c - b))
    e2 = norm(np.cross(plane_normal, e1))
    return c + bond_length * (
        math.cos(bond_angle) * e1
        + math.sin(bond_angle) * (math.cos(phi) * e2 + math.sin(phi) * plane_normal)
    )

def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    value = np.dot(norm(a - b), norm(c - b))
    return math.degrees(math.acos(float(np.clip(value, -1.0, 1.0))))

def dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    b0 = a - b
    b1 = c - b
    b2 = d - c
    b1n = norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))
