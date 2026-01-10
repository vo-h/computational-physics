from scipy.integrate import tplquad
import math
from pydantic import BaseModel, Field
from source.hf.sto_ng import STOGOrbital, SijIntegrator, TijIntegrator
from source.hf import L
import numpy as np
from tqdm import tqdm
from typing import Literal, Callable
from source.hf.atom import Atom
import numpy as np
from functools import partial

overlap_integrator = SijIntegrator()
kinetic_integrator = TijIntegrator()

class Molecule(BaseModel):
    """Molecule class, which consists of multiple atoms."""
    atoms: list[Atom]
    cij: list[float] = Field(default_factory=list, description="The coefficients to be optimized during the SCF procedure.")

    @property
    def orbitals(self) -> list[STOGOrbital]:
        """List of all atomic orbitals in the molecule."""
        return [orb for atom in self.atoms for orb in atom.orbitals]

    @property
    def cc(self) -> np.ndarray:
        """The contracted coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.cc for orb in self.orbitals])
    
    @property
    def alpha(self) -> np.ndarray:
        """The exponent coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.alpha for orb in self.orbitals])
    
    def __call__(self, x: float, y: float, z: float) -> np.ndarray:
        """Evaluate the orbitals of the atom at a given point (x, y, z)."""
        return np.array([orb(x, y, z) for orb in self.orbitals])
    
    def eval_ao(self, coords: np.ndarray) -> np.ndarray:
        """Evaluate the orbitals of the molecule at a given set of coordinates."""
        return np.array([[orb(x, y, z) for orb in self.orbitals] for x, y, z in coords])
    
    @property
    def S(self) -> np.ndarray:
        """
        The overlap integral between all pairs of orbitals in the molecule.
        https://content.wolfram.com/sites/19/2012/02/Ho.pdf
        """
        S = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                if i == j:
                    S[i, j] = 1.0  # The overlap of an orbital with itself is 1
                else:
                    S[i, j] = S[j, i] = self._Sij(i, j)

        return S

    @property
    def T(self) -> np.ndarray:
        """The kinetic energy integral between all pairs of orbitals in the molecule."""
        T = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                Tij = self._Tij(i, j)
                T[i, j] = self._Tij(i, j)
                if i != j:
                    T[j, i] = Tij  # The kinetic energy integral is symmetric
        return T

    @property
    def V(self) -> np.ndarray:
        """The electron-nuclear attraction integral between all pairs of orbitals in the molecule."""
        V = np.zeros((len(self.orbitals), len(self.orbitals)))
        pbar = tqdm(total=len(self.orbitals)**2)
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                func = lambda z, y, x: self.orbitals[i](x=x, y=y, z=z) * self.orbitals[j](x=x, y=y, z=z) * sum(-atom.atom_charge / math.sqrt((x - atom.coords[0])**2 + (y - atom.coords[1])**2 + (z - atom.coords[2])**2) for atom in self.atoms)
                V[i, j] = V[j, i] = tplquad(func, -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)[0]
                pbar.update(2)
        return V

    @property
    def H(self) -> np.ndarray:
        """The core Hamiltonian matrix, which is the sum of the kinetic energy and electron-nuclear attraction integrals."""
        return self.T + self.V

    @property
    def J(self) -> np.ndarray:
        """The Coulomb matrix, which is calculated from the density matrix and the two-electron integrals."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the Coulomb matrix from the density matrix and the two-electron integrals.
        return np.zeros((len(self.orbitals), len(self.orbitals)))

    @property
    def K(self) -> np.ndarray:
        """The exchange matrix, which is calculated from the density matrix and the two-electron integrals."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the exchange matrix from the density matrix and the two-electron integrals.
        return np.zeros((len(self.orbitals), len(self.orbitals)))

    @property
    def P(self) -> np.ndarray:
        """The density matrix, which is calculated from the coefficients of the molecular orbitals."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the density matrix from the molecular orbital coefficients.
        return np.zeros((len(self.orbitals), len(self.orbitals)))

    @property
    def F(self) -> np.ndarray:
        """The Fock matrix, which is the sum of the core Hamiltonian and the electron-electron repulsion integrals."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the Fock matrix from the core Hamiltonian and the electron-electron repulsion integrals.
        return self.H + np.zeros((len(self.orbitals), len(self.orbitals)))

    def _Sij(self, i: int, j: int, resolver: Literal['Python', 'C'] = 'Python') -> float:
        """Calculate the overlap integral between two orbitals."""
        if resolver == "Python":
            return overlap_integrator.Sij(self.orbitals[i], self.orbitals[j])
        return 0.0  # Placeholder return value for C implementation, which should be replaced with the actual calculation of the overlap integral using the C library.

    def _Tij(self, i: int, j: int, resolver: Literal['Python', 'C'] = 'Python') -> float:
        """Calculate the kinetic energy integral between two orbitals."""
        if resolver == "Python":
            return kinetic_integrator.Tij(self.orbitals[i], self.orbitals[j])
        return 0.0  # Placeholder return value for C implementation, which should be replaced with the actual calculation of the kinetic energy integral using the C library.

    def to_pyscf_coords(self) -> str:   
        """Convert the molecule's atomic coordinates to a format compatible with PySCF."""
        return ";".join(f"{atom.atom} {atom.coords[0]} {atom.coords[1]} {atom.coords[2]}" for atom in self.atoms)
    