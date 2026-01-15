from xml.dom.minidom import Element
from scipy.integrate import tplquad
import math
from pydantic import BaseModel, ConfigDict, Field, model_validator
from hf.basis_stog import STOGOrbital, STOGIntegrator
from hf import L
import numpy as np
from tqdm import tqdm
from hf.atom import Atom
import numpy as np

class Molecule(BaseModel):
    """Molecule class, which consists of multiple atoms."""
    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)
    atoms: list[Atom]
    cij: list[float] = Field(default_factory=list, description="The coefficients to be optimized during the SCF procedure.")

    @model_validator(mode="after")
    def validate_basis_consistency(self):
        bases = set(atom.basis for atom in self.atoms)
        if len(bases) > 1:
            raise ValueError(f"All atoms in the molecule must have the same basis set. Found the following basis sets: {', '.join(bases)}.")
        return self

    @property
    def orbitals(self) -> list[STOGOrbital]:
        """List of all atomic orbitals in the molecule."""
        return [orb for atom in self.atoms for orb in atom.orbitals]

    @property
    def cc(self) -> np.ndarray:
        """The contracted coefficients of the GTOs in the STO-nG expansion as a numpy array."""
        return np.stack([orb.cc for orb in self.orbitals])
    
    @property
    def alpha(self) -> np.ndarray:
        """The exponent coefficients of the GTOs in the STO-nG expansion as a numpy array."""
        return np.stack([orb.alpha for orb in self.orbitals])
    
    @property
    def basis(self) -> str:
        """Return a list of all the STOGOrbitals in the molecule."""
        return self.atoms[0].basis
    
    @property
    def integrator(self) -> STOGIntegrator:
        """Return an instance of the STOGIntegrator class, which can be used to calculate the matrices."""
        if self.basis.startswith("STO-"):
            return STOGIntegrator()
    
    @property
    def S(self) -> np.ndarray:
        """The overlap integral between all pairs of orbitals in the molecule."""
        S = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                S[i, j] = S[j, i] = self.integrator.Sij(self.orbitals[i], self.orbitals[j])
        return S

    @property
    def T(self) -> np.ndarray:
        """The kinetic energy integral between all pairs of orbitals in the molecule."""
        T = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                T[i, j] = T[j, i] = self.integrator.Tij(self.orbitals[i], self.orbitals[j])
        return T

    @property
    def V(self) -> np.ndarray:
        """The electron-nuclear attraction integral between all pairs of orbitals in the molecule."""
        pbar = tqdm(total=len(self.orbitals)*len(self.orbitals), desc="Calculating V matrix")
        V = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                term = sum(-atom.Z * self.integrator.VijR(self.orbitals[i], self.orbitals[j], atom.coords) for atom in self.atoms)
                V[i, j] = V[j, i] = term
                pbar.update(1 if i == j else 2)
        pbar.close()
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

    def to_pyscf_coords(self) -> str:   
        """Convert the molecule's atomic coordinates to a format compatible with PySCF."""
        return ";".join(f"{atom.atom} {atom.coords[0]} {atom.coords[1]} {atom.coords[2]}" for atom in self.atoms)
    
    def __call__(self, x: float, y: float, z: float) -> np.ndarray:
        """Evaluate the orbitals of the atom at a given point (x, y, z)."""
        return np.array([orb(x, y, z) for orb in self.orbitals])
    
    def HF(self):
        """Perform the Hartree-Fock self-consistent field (SCF) procedure to optimize the coefficients of the molecular orbitals."""
        # This is a placeholder implementation and should be replaced with the actual implementation of the Hartree-Fock SCF procedure.
        pass

    def optimize(self):
        """Optimize the geometry of the molecule."""
        # This is a placeholder implementation and should be replaced with the actual implementation of an optimization algorithm to optimize the geometry of the molecule.
        pass
