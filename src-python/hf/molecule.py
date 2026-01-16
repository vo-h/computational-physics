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
from functools import cached_property

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
    
    @cached_property
    def S(self) -> np.ndarray:
        """The overlap integral between all pairs of orbitals in the molecule."""
        S = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                S[i, j] = S[j, i] = self.integrator.Sij(self.orbitals[i], self.orbitals[j])
        return S
    
    @cached_property
    def S12(self) -> np.ndarray:
        """The symmetric orthogonalization matrix"""
        eigv, L = np.linalg.eig(self.S)
        LAMBDA12 = np.diag(1/np.sqrt(eigv))
        res = L @ LAMBDA12 @ np.transpose(L)
        return np.where(np.isclose(res, 0.0), 0.0, res)

    @cached_property
    def T(self) -> np.ndarray:
        """The kinetic energy integral between all pairs of orbitals in the molecule."""
        T = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                T[i, j] = T[j, i] = self.integrator.Tij(self.orbitals[i], self.orbitals[j])
        return T

    @cached_property
    def Vne(self) -> np.ndarray:
        """The electron-nuclear attraction integral between all pairs of orbitals in the molecule."""
        pbar = tqdm(total=len(self.orbitals)*len(self.orbitals), desc="Calculating V matrix")
        V = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                V[i, j] = V[j, i] = sum(-atom.Z * self.integrator.VijR(self.orbitals[i], self.orbitals[j], atom.coords) for atom in self.atoms)
                pbar.update(1 if i == j else 2)
        return V

    @cached_property
    def H(self) -> np.ndarray:
        """The core Hamiltonian matrix, which is the sum of the kinetic energy and electron-nuclear attraction integrals."""
        return self.T + self.Vne

    @cached_property
    def Vee(self) -> np.ndarray:
        """The electron-electron repulsion integrals (ij|kl) with 8-fold permutational symmetry."""
        V = np.zeros((len(self.orbitals), len(self.orbitals), len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(len(self.orbitals)):
                for k in range(len(self.orbitals)):
                    for l in range(len(self.orbitals)):
                        if (i*len(self.orbitals) + j) >= (k*len(self.orbitals) + l):
                            term = self.integrator.Vijkl(self.orbitals[i], self.orbitals[j], self.orbitals[k], self.orbitals[l])
                            # 8-fold symmetry: (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (kl|ji) = (lk|ij) = (lk|ji)
                            V[i,j,k,l] = term
                            V[j,i,k,l] = term
                            V[i,j,l,k] = term
                            V[j,i,l,k] = term
                            V[k,l,i,j] = term
                            V[k,l,j,i] = term
                            V[l,k,i,j] = term
                            V[l,k,j,i] = term
        return V


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

    @cached_property
    def D0(self) -> float:
        """The initial guess for the total electronic energy of the molecule, which is calculated from the core Hamiltonian and the density matrix."""
        # F0 = np.transpose(self.S12) @ self.H @ self.S12
        # eigv, C0 = np.linalg.eig(F0)
        # C0 = self.S12 @ C0
        # return C0

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
