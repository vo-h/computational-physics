from scipy.integrate import tplquad
import math
from pydantic import BaseModel
from functools import cached_property
from source.hartree_fock.basis_sets import STO3GOrbital, Orbital
from source.hartree_fock import L, STO3G
import numpy as np
from tqdm import tqdm
from typing import Literal

class Atom(BaseModel):
    """STO-3G for an atom, which consists of multiple orbitals."""
    atom: str
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)
    basis_set: Literal["STO-3G"] = "STO-3G"

    @cached_property
    def orbitals(self) -> list[Orbital]:
        orbitals = []
        occupations = {"S": 1, "P": 3, "D": 5, "F": 7}

        for orbital in STO3G[self.atom]['cgfs']:
            occ = occupations[orbital['type']]
            for ind in range(occ):
                if self.basis_set == "STO-3G":
                    orbitals.append(STO3GOrbital(
                        type=orbital['type'],
                        atom=self.atom,
                        coordinates=self.coordinates,
                        params1=(orbital['gtos'][0]['coeff'], orbital['gtos'][0]['alpha'], L[orbital['type']][ind][0], L[orbital['type']][ind][1], L[orbital['type']][ind][2]),
                        params2=(orbital['gtos'][1]['coeff'], orbital['gtos'][1]['alpha'], L[orbital['type']][ind][0], L[orbital['type']][ind][1], L[orbital['type']][ind][2]),
                        params3=(orbital['gtos'][2]['coeff'], orbital['gtos'][2]['alpha'], L[orbital['type']][ind][0], L[orbital['type']][ind][1], L[orbital['type']][ind][2]),
                    ))
        return orbitals
    
    @cached_property
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.prime_coefficients for orb in self.orbitals])

    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.orbital_coefficients for orb in self.orbitals])

    @cached_property
    def normalization_constants(self) -> np.ndarray:
        """Calculate the normalization constant for the STO-3G orbital."""
        return np.stack([orb.normalization_constants for orb in self.orbitals])

class Molecule(BaseModel):
    """STO-3G for a molecule, which consists of multiple atoms."""
    atoms: list[Atom]

    @cached_property
    def orbitals(self) -> list[Orbital]:
        return [orb for atom in self.atoms for orb in atom.orbitals]
    
    @cached_property
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.prime_coefficients for orb in self.orbitals])

    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.orbital_coefficients for orb in self.orbitals])

    @cached_property
    def normalization_constants(self) -> np.ndarray:
        """Calculate the normalization constant for the STO-3G orbital."""
        return np.stack([orb.normalization_constants for orb in self.orbitals])

    def overlap(self) -> float:
        """
        Calculate the overlap integral between all pairs of orbitals in the molecule.
        https://content.wolfram.com/sites/19/2012/02/Ho.pdf
        """
        S = np.zeros((len(self.orbitals), len(self.orbitals)))
        pbar = tqdm(total=len(self.orbitals)**2)
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                if i == j:
                    S[i, j] = 1.0  # The overlap of an orbital with itself is 1
                    pbar.update(1)
                else:
                    func = lambda z, y, x: self.orbitals[i](x=x, y=y, z=z) * self.orbitals[j](x=x, y=y, z=z)
                    s = tplquad(func, -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)[0]
                    S[i, j] = S[j, i] = 0.0 if math.isclose(s, 0.0, abs_tol=1e-10) else s  # Set very small overlaps to zero to avoid numerical issues
                    pbar.update(2)
        return S
    
    def kinetic_energy(self) -> float:
        """Calculate the kinetic energy integral between all pairs of orbitals in the molecule."""
        T = np.zeros((len(self.orbitals), len(self.orbitals)))
        pbar = tqdm(total=len(self.orbitals)**2)
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                func = lambda z, y, x: self.orbitals[i](x=x, y=y, z=z) * self.orbitals[j](x=x, y=y, z=z) * (-0.5 * (self.orbitals[j].GTO1.alpha * 3 - 2 * self.orbitals[j].GTO1.alpha**2 * ((x - self.orbitals[j].GTO1.coordinates[0])**2 + (y - self.orbitals[j].GTO1.coordinates[1])**2 + (z - self.orbitals[j].GTO1.coordinates[2])**2)))
                T[i, j] = T[j, i] = tplquad(func, -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)[0]
                pbar.update(2)
        return T
    
    



    



    



