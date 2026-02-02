from functools import cached_property
import math
from pydantic import BaseModel, ConfigDict, Field, model_validator
from hf.basis_stog import STOGOrbital, STOGIntegrator, STOGDerivator
from hf import L
import numpy as np
from tqdm import tqdm
from hf.atom import Atom
import numpy as np
import logging
from numpy.typing import NDArray
from typing import Self, Optional

logger = logging.getLogger(__name__)

class Molecule(BaseModel):
    """Molecule class, which consists of multiple atoms."""
    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)
    atoms: list[Atom]
    C: Optional[NDArray[np.float64]] = Field(default=None, description="The coefficients to be optimized during the SCF procedure.")
    F: Optional[NDArray[np.float64]] = Field(default=None, description="The Fock matrix to be optimized during the SCF procedure.")
    E: Optional[float] = Field(default=None, description="The electronic energy of the molecule, which is optimized during the SCF procedure.")

    @model_validator(mode="after")
    def validate_basis_consistency(self):
        bases = set(atom.basis for atom in self.atoms)
        if len(bases) > 1:
            raise ValueError(f"All atoms in the molecule must have the same basis set. Found the following basis sets: {', '.join(bases)}.")
        return self
    
    @model_validator(mode="after")
    def check_electron_count(self):
        total_electrons = sum(atom.Z for atom in self.atoms)
        if total_electrons % 2 != 0:
            raise ValueError(f"The molecule has an odd number of electrons ({total_electrons}), which is not supported by the current implementation. Please use a molecule with an even number of electrons.")
        return self

    @classmethod
    def from_file(cls, file_path: str) -> Self:
        """Create a Molecule object from an input file. The file should contain the atomic symbols and coordinates in the following format:
        
        C   0.000000    0.000000    0.000000    6   STO-3G
        H   0.000000    0.000000    1.089000    1   STO-3G
        H   1.026719    0.000000   -0.363000    1   STO-3G
        H  -0.513360   -0.889165   -0.363000    1   STO-3G
        H  -0.513360    0.889165   -0.363000    1   STO-3G
        
        Each line corresponds to an atom, with the atomic symbol followed by the x, y, z coordinates in Angstroms and a basis set.
        """
        atoms = []
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) != 6:
                    raise ValueError(f"Invalid line in input file: {line.strip()}. Each line must contain an atomic symbol followed by three coordinates and a basis set.")
                atom_symbol = parts[0]
                coords = tuple(float(coord) for coord in parts[1:4])
                basis_set = parts[5]
                atoms.append(Atom(atom=atom_symbol, coords=coords, basis_set=basis_set))
        return cls(atoms=atoms)


    @property
    def orbitals(self) -> list[STOGOrbital]:
        """List of all atomic orbitals in the molecule."""
        return [orb for atom in self.atoms for orb in atom.orbitals]
    
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
    def derivator(self) -> STOGDerivator:
        """Return an instance of the STOGDerivator class, which can be used to calculate the derivative matrices."""
        if self.basis.startswith("STO-"):
            return STOGDerivator()
    
    @cached_property
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
    def Vne(self) -> np.ndarray:
        """The electron-nuclear attraction integral between all pairs of orbitals in the molecule."""
        V = np.zeros((len(self.orbitals), len(self.orbitals)))
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                V[i, j] = V[j, i] = sum(-atom.Z * self.integrator.VijR(self.orbitals[i], self.orbitals[j], atom.coords) for atom in self.atoms)
        return V

    @property
    def H(self) -> np.ndarray:
        """The core Hamiltonian matrix, which is the sum of the kinetic energy and electron-nuclear attraction integrals."""
        return self.T + self.Vne

    @property
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
    def Vnn(self) -> float:
        """The nuclear-nuclear repulsion energy of the molecule."""
        Vnn = 0.0
        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                R = np.linalg.norm(np.array(self.atoms[i].coords) - np.array(self.atoms[j].coords))
                Vnn += self.atoms[i].Z * self.atoms[j].Z / R
        return Vnn

    @property
    def S12(self) -> np.ndarray:
        """The symmetric orthogonalization matrix"""
        eigv, L = np.linalg.eig(self.S)
        LAMBDA12 = np.diag(1/np.sqrt(eigv))
        res = L @ LAMBDA12 @ L.T
        return np.where(np.isclose(res, 0.0), 0.0, res)

    @property
    def F0(self) -> np.ndarray:
        """The initial guess for the Fock matrix, which is the core Hamiltonian."""
        F0 = self.S12.T @ self.H @ self.S12
        return np.where(np.isclose(F0, 0.0), 0.0, F0)

    @property
    def C0(self) -> np.ndarray:
        """The initial guess for the coefficients of the molecular orbitals."""
        eigv, C0 = np.linalg.eig(self.F0)
        idx = np.argsort(eigv) # Order the eigenvalues and eigenvectors in ascending order of eigenvalues
        eigv, C0 = eigv[idx], C0[:, idx]
        C0 = self.S12 @ C0 # Transform the eigenvectors back to the original basis
        return np.where(np.isclose(C0, 0.0), 0.0, C0)

    @property
    def D0(self) -> np.ndarray:
        """The initial guess for the density matrix."""
        occ = math.ceil(sum(atom.Z for atom in self.atoms) / 2)
        return self.C0[:, :occ] @ self.C0[:, :occ].T
    
    @property
    def E0(self) -> float:
        """The initial guess for the electronic energy of the molecule.
        
        For the initial SCF electronic energy: E0 = 2 * Σ_μν D[μ,ν] * H[μ,ν]
        The factor of 2 accounts for spin (each spatial orbital holds 2 electrons).
        """
        return 2 * np.sum(self.D0 * self.H)

    # @property
    # def gradients(self) -> np.ndarray:
    #     """Calculate the nuclear gradients of the molecule using analytical derivatives."""
    #     return np.zeros((len(self.atoms), 3))
    
    @property
    def gradients(self) -> np.ndarray:
        """Calculate the nuclear gradients of the molecule using analytical derivatives.

        Implements the standard RHF analytic gradient (Pulay + 1e/2e + nuclear repulsion) in the AO basis:

            dE/dR_A = sum_{μν} P_{μν} dH_{μν}/dR_A
                    + 1/2 sum_{μνλσ} P_{μν} P_{λσ} d(μν|λσ)/dR_A
                    - sum_{μν} W_{μν} dS_{μν}/dR_A
                    + dV_nn/dR_A

        where P is the spin-summed density (P = 2D for closed-shell RHF) and
        W is the spin-summed energy-weighted density (W = 2 C_occ ε_occ C_occ^T).
        """
        # Ensure we have a converged SCF solution
        if getattr(self, "C", None) is None or getattr(self, "F", None) is None or self.C is None or self.F is None or self.C.size == 0:
            self.CHF()

        nbf = len(self.orbitals)
        nat = len(self.atoms)
        occ = math.ceil(sum(atom.Z for atom in self.atoms) / 2)

        # Recompute MO energies ε (needed for W) from the converged Fock matrix
        S12 = self.S12
        eigv, C_orth = np.linalg.eig(S12.T @ self.F @ S12)
        idx = np.argsort(eigv)
        eigv, C_orth = eigv[idx], C_orth[:, idx]
        C = S12 @ C_orth

        # Density matrices (closed-shell RHF, same convention as CHF())
        D = C[:, :occ] @ C[:, :occ].T   # un-doubled
        P = 2.0 * D                     # spin-summed density

        # Energy-weighted density (spin-summed)
        W = 2.0 * (C[:, :occ] @ np.diag(eigv[:occ]) @ C[:, :occ].T)

        deriv = self.derivator
        dirs = ("x", "y", "z")

        grad = np.zeros((nat, 3))

        # Atom-label string used by orbitals (orb.atom) must match what we pass into derivator
        atom_labels = [atom.atom for atom in self.atoms]
        pbar = tqdm(total=nat*3, desc="Computing Gradients", unit="term")
        for A, atomA in enumerate(self.atoms):
            RA = np.array(atomA.coords, dtype=float)
            ZA = float(atomA.Z)

            # Nuclear-nuclear repulsion gradient for atom A
            g_nn = np.zeros(3)
            for B, atomB in enumerate(self.atoms):
                if B == A:
                    continue
                RB = np.array(atomB.coords, dtype=float)
                ZB = float(atomB.Z)
                rvec = RA - RB
                r = np.linalg.norm(rvec)
                if r > 0:
                    g_nn -= ZA * ZB * rvec / (r**3)  # Note: negative sign!

            for a, d in enumerate(dirs):
                # Build AO-derivative matrices dS, dT, dVne (explicit nucleus term)
                dS = np.zeros((nbf, nbf))
                dT = np.zeros((nbf, nbf))
                dVne = np.zeros((nbf, nbf))

                for i in range(nbf):
                    oi = self.orbitals[i]
                    for j in range(i, nbf):
                        oj = self.orbitals[j]

                        # AO-center derivatives (angular-momentum shift identity)
                        dS_ij = deriv.dSij_or_dTij(oi, oj, atom_labels[A], d, matrix="S")
                        dT_ij = deriv.dSij_or_dTij(oi, oj, atom_labels[A], d, matrix="T")

                        # Explicit derivative of e–n attraction wrt nucleus A:
                        # Testing negative sign for field integral
                        dVne_ij = -ZA * deriv.dVijR(oi, oj, atomA.coords, d)

                        dS[i, j] = dS[j, i] = dS_ij
                        dT[i, j] = dT[j, i] = dT_ij
                        dVne[i, j] = dVne[j, i] = dVne_ij

                dH = dT + dVne

                # One-electron term: 2 * Tr[D dH] (factor of 2 for closed-shell)
                g1 = 2.0 * np.sum(D * dH)

                # Pulay term: -Tr[W dS]
                g_pulay = -np.sum(W * dS)

                # Two-electron term: Sum D_ij D_kl [2 d(ij|kl) - d(ik|jl)]
                # Using un-doubled density D for proper formula
                g2 = 0.0
                for i in range(nbf):
                    for j in range(nbf):
                        Dij = D[i, j]
                        if abs(Dij) < 1e-14:
                            continue
                        for k in range(nbf):
                            for l in range(nbf):
                                Dkl = D[k, l]
                                if abs(Dkl) < 1e-14:
                                    continue
                                
                                # Coulomb derivative: 2 * d(ij|kl)
                                dERI_coulomb = deriv.dVijkl(
                                    self.orbitals[i], self.orbitals[j],
                                    self.orbitals[k], self.orbitals[l],
                                    atom_labels[A], d
                                )
                                
                                # Exchange derivative: -d(ik|jl)
                                dERI_exchange = deriv.dVijkl(
                                    self.orbitals[i], self.orbitals[k],
                                    self.orbitals[j], self.orbitals[l],
                                    atom_labels[A], d
                                )
                                
                                g2 += Dij * Dkl * (2.0 * dERI_coulomb - dERI_exchange)

                grad[A, a] = g1 + g2 + g_pulay + g_nn[a]
                pbar.update(1)
        pbar.close()

        return np.where(np.isclose(grad, 0.0), 0.0, grad)


    @property
    def hessian(self) -> np.ndarray:    
        """Calculate the nuclear Hessian of the molecule."""
        return np.zeros((len(self.atoms), 3, 3))

    def to_pyscf_coords(self) -> str:   
        """Convert the molecule's atomic coordinates to a format compatible with PySCF."""
        return ";".join(f"{atom.atom} {atom.coords[0]} {atom.coords[1]} {atom.coords[2]}" for atom in self.atoms)

    def __call__(self, x: float, y: float, z: float, mo: int) -> float:
        """Evaluate a molecular orbital of the atom at a given point (x, y, z)."""
        if self.C is None:
            raise ValueError(
                "SCF calculation has not been performed yet. "
                "Please run SCF before evaluating molecular orbitals."
            )
        
        val = 0.0
        for i, orbital in enumerate(self.orbitals):
            val += self.C[i, mo] * orbital(x, y, z)
        return val
    
    def CHF(self, delta: float = 1e-12, max_iter: int = 100):
        """Perform the Hartree-Fock self-consistent field (SCF) procedure to optimize the coefficients of the molecular orbitals."""
        
        occ = math.ceil(sum(atom.Z for atom in self.atoms) / 2)
        D = self.D0
        E = self.E0
        H = self.H
        Vee = self.Vee
        S12 = self.S12

        for i in tqdm(range(max_iter), desc="SCF Iteration"):
            if i == 0:
                logger.info(f"Iteration 00: Energy = {E:.12f} Hartree")
                continue
            F = H + np.einsum('ls,uvls->uv', D, 2*Vee) - np.einsum('ls,ulvs->uv', D, Vee) # Compute new Fock matrix in AO basis using tensor contraction
            eigv, C = np.linalg.eig(S12.T @ F @ S12) # Orthogonalize & Diagonalize Fock matrix
            idx = np.argsort(eigv) # Sort eigenvalues and eigenvectors
            eigv, C = eigv[idx], C[:, idx]
            C = S12 @ C # Transform eigenvectors back to original basis
            D = C[:, :occ] @ C[:, :occ].T # Build new density matrix
            E_new = np.sum(D * (H + F)) # Compute new energy
            if abs(E_new - E) < delta:
                logger.info(f"SCF converged in {i} iterations with energy: {E_new:.6f} Hartree")
                break
            logger.info(f"Iteration {i}: Energy = {E_new:.6f} Hartree, ΔE = {E_new - E:.6e} Hartree")
            E = E_new
        self.C = C
        self.F = F
        self.E = E_new

    def optimize(self):
        """Optimize the geometry of the molecule."""
        pass
    