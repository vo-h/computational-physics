from scipy.integrate import tplquad
import math
from pydantic import BaseModel
from functools import cached_property
from source.hartree_fock.basis_sets import STO3G
from source.hartree_fock import L
import numpy as np
from tqdm import tqdm


class GTO(BaseModel):
    """Gaussian Type Orbital (GTO) basis function."""
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)

    # Exponent and coefficient for the GTO.
    c: float = 1.0
    alpha: float

    # Angular momentum quantum numbers for the GTO (s, p, d, etc.)
    nx: int = 0
    ny: int = 0
    nz: int = 0

    @cached_property
    def N(self) -> float:
        """Normalization constant for the GTO. https://en.wikipedia.org/wiki/Gaussian_orbital#Cartesian_coordinates"""
        prefactor = (2 * self.alpha / math.pi) ** (3 / 4)
        numerator = (8 * self.alpha) ** (self.nx + self.ny + self.nz) * math.factorial(self.nx) * math.factorial(self.ny) * math.factorial(self.nz)
        denominator = math.factorial(2*self.nx) * math.factorial(2*self.ny) * math.factorial(2*self.nz)
        return prefactor * math.sqrt(numerator/ denominator)

    def __call__(self, z: float, y: float, x: float) -> float:
        """Evaluate the GTO at a given point (x, y, z)."""
        r_squared = (x - self.coordinates[0]) ** 2 + (y - self.coordinates[1]) ** 2 + (z - self.coordinates[2]) ** 2
        return self.N * (x - self.coordinates[0]) ** self.nx * (y - self.coordinates[1]) ** self.ny * (z - self.coordinates[2]) ** self.nz * math.exp(-self.alpha * r_squared)


class STO3GOrbital(BaseModel):
    """STO-3G for an orbital (1s, 2s, 2p, etc.), which is a linear combination of 3 GTOs."""
    type: str
    atom: str | None = None
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Exponent and coefficient for GTOs in the STO-3G expansion
    params1: tuple[float, float, int, int, int]  # (c1, alpha1, nx1, ny1, nz1)
    params2: tuple[float, float, int, int, int]  # (c2, alpha2, nx2, ny2, nz2)
    params3: tuple[float, float, int, int, int]  # (c3, alpha3, nx3, ny3, nz3)

    @cached_property
    def GTO1(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params1[1], nx=self.params1[2], ny=self.params1[3], nz=self.params1[4])
    
    @cached_property
    def GTO2(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params2[1], nx=self.params2[2], ny=self.params2[3], nz=self.params2[4])
    
    @cached_property
    def GTO3(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params3[1], nx=self.params3[2], ny=self.params3[3], nz=self.params3[4])

    @cached_property
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[0], self.params2[0], self.params3[0]])

    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[1], self.params2[1], self.params3[1]])

    def __call__(self, z: float, y: float, x: float) -> float:
        """Evaluate the STO-3G orbital at a given point (x, y, z)."""
        gto1 = self.GTO1(x, y, z)
        gto2 = self.GTO2(x, y, z)
        gto3 = self.GTO3(x, y, z)
        return self.params1[0] * gto1 + self.params2[0] * gto2 + self.params3[0] * gto3
    
class STO3GAtom(BaseModel):
    """STO-3G for an atom, which consists of multiple orbitals."""
    atom: str
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)

    @cached_property
    def orbitals(self) -> list[STO3GOrbital]:
        orbitals = []
        occupations = {"S": 1, "P": 3, "D": 5, "F": 7}

        for orbital in STO3G[self.atom]['cgfs']:
            occ = occupations[orbital['type']]
            for ind in range(occ):
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

class STO3Molecule(BaseModel):
    """STO-3G for a molecule, which consists of multiple atoms."""
    atoms: list[STO3GAtom]

    @cached_property
    def orbitals(self) -> list[STO3GOrbital]:
        return [orb for atom in self.atoms for orb in atom.orbitals]
    
    @cached_property
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.prime_coefficients for orb in self.orbitals])

    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.orbital_coefficients for orb in self.orbitals])

    def overlap(self) -> float:
        """Calculate the overlap integral between all pairs of orbitals in the molecule."""
        S = -np.ones((len(self.orbitals), len(self.orbitals)))
        pbar = tqdm(total=len(self.orbitals)**2)
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                if i == j:
                    S[i, j] = 1.0  # The overlap of an orbital with itself is 1
                else:
                    func = lambda z, y, x: self.orbitals[i](z, y, x) * self.orbitals[j](z, y, x)
                    S[i, j] = S[j, i] = tplquad(func, -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)[0]
                pbar.update(2)
        return S
    
    def kinetic_energy(self) -> float:
        """Calculate the kinetic energy integral between all pairs of orbitals in the molecule."""
        T = -np.ones((len(self.orbitals), len(self.orbitals)))
        pbar = tqdm(total=len(self.orbitals)**2)
        for i in range(len(self.orbitals)):
            for j in range(i, len(self.orbitals)):
                func = lambda z, y, x: self.orbitals[i](z, y, x) * self.orbitals[j](z, y, x) * (-0.5 * (self.orbitals[j].GTO1.alpha * 3 - 2 * self.orbitals[j].GTO1.alpha**2 * ((x - self.orbitals[j].GTO1.coordinates[0])**2 + (y - self.orbitals[j].GTO1.coordinates[1])**2 + (z - self.orbitals[j].GTO1.coordinates[2])**2)))
                T[i, j] = T[j, i] = tplquad(func, -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)[0]
                pbar.update(2)
        return T
    
    



    



    



