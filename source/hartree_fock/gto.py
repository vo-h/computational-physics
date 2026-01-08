import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field, computed_field
from functools import cached_property
import numpy as np
from abc import ABC, abstractmethod
from scipy.integrate import tplquad
from scipy import special
from typing import Callable

class GTO(BaseModel):
    """Gaussian Type Orbital (GTO) basis function."""
    coords: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)

    # Exponent and angular momentum quantum numbers for the GTO (s, p, d, etc.)
    alpha: float
    nx: int = 0
    ny: int = 0
    nz: int = 0

    @cached_property
    def l(self) -> int:
        """Calculate the total angular momentum quantum number."""
        return self.nx + self.ny + self.nz

    @cached_property
    def N_cartesian(self) -> float:
        """Normalization constant for the cartisian GTO. https://en.wikipedia.org/wiki/Gaussian_orbital#Cartesian_coordinates"""
        prefactor = (2 * self.alpha / math.pi) ** (3 / 4)
        numerator = (8 * self.alpha) ** (self.nx + self.ny + self.nz) * math.factorial(self.nx) * math.factorial(self.ny) * math.factorial(self.nz)
        denominator = math.factorial(2*self.nx) * math.factorial(2*self.ny) * math.factorial(2*self.nz)
        return prefactor * math.sqrt(numerator/ denominator)
    
    @cached_property
    def N_spherical(self) -> float:
        """Calculate the normalization constant for the radial/spherical GTO. https://pyscf.org/pyscf_api_docs/pyscf.gto.html#pyscf.gto.mole.MoleBase.gto_norm"""
        num_1 = 2**(2*self.l+3)
        num_2 = math.factorial(self.l+1)
        num_3 = (2*self.alpha)**(self.l+1.5)
        denom = math.sqrt(math.pi) * math.factorial(2*self.l+2)
        return math.sqrt(num_1 * num_2 * num_3 / denom)

    @cached_property
    def Y_spherical(self) -> Callable[[float, float, float], float]:
        """Calculate the angular part of the GTO. https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics"""

        # For s orbitals
        if self.l == 0:
            return lambda x, y, z: 0.5 * math.sqrt(1 / math.pi)
        
        # For p orbitals
        if self.l == 1:
            prefactor = lambda x, y, z: math.sqrt(0.75 / math.pi) / self.dist(x, y, z)
            if self.nx == 1:
                return lambda x, y, z: prefactor(x, y, z) * (x - self.coords[0])
            elif self.ny == 1:
                return lambda x, y, z: prefactor(x, y, z) * (y - self.coords[1])
            elif self.nz == 1:
                return lambda x, y, z: prefactor(x, y, z) * (z - self.coords[2])
        
        # For d orbitals
        if self.l == 2:
            prefactor = lambda x, y, z: math.sqrt(15 / (4 * math.pi)) / self.r_squared(x, y, z)
            if self.nx == 2 or self.ny == 2: # d_x2-y2
                return lambda x, y, z: 0.5 * prefactor(x, y, z) * ((x - self.coords[0])**2 - (y - self.coords[1])**2)
            elif self.nz == 2: # d_z2
                return lambda x, y, z: 0.5 * prefactor(x, y, z) * (3 * (z - self.coords[2])**2 - self.r_squared(x, y, z))
            elif self.nx == 1 and self.ny == 1: # d_xy
                return lambda x, y, z: prefactor(x, y, z) * (x - self.coords[0]) * (y - self.coords[1])
            elif self.nx == 1 and self.nz == 1: # d_xz
                return lambda x, y, z: prefactor(x, y, z) * (x - self.coords[0]) * (z - self.coords[2])
            elif self.ny == 1 and self.nz == 1: # d_yz
                return lambda x, y, z: prefactor(x, y, z) * (y - self.coords[1]) * (z - self.coords[2])

    @cached_property
    def Y_cartesian(self) -> Callable[[float, float, float], float]:
        """Calculate the angular part of the GTO in cartesian coordinates."""
        return lambda x, y, z: (x - self.coords[0]) ** self.nx * (y - self.coords[1]) ** self.ny * (z - self.coords[2]) ** self.nz

    def r_squared(self, x: float, y: float, z: float) -> float:
        """Calculate the squared distance from the center of the GTO."""
        dist_x = x - self.coords[0]
        dist_y = y - self.coords[1]
        dist_z = z - self.coords[2]
        return dist_x ** 2 + dist_y ** 2 + dist_z ** 2

    def dist(self, x: float, y: float, z: float) -> float:
        """Calculate the distance from the center of the GTO."""
        return math.sqrt(self.r_squared(x, y, z))

    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the GTO at a given point (x, y, z)."""
        R = math.exp(-self.alpha * self.r_squared(x, y, z))
        return self.N_cartesian * R * self.Y_cartesian(x, y, z) # Use cartesian coordinates


class Orbital(BaseModel, ABC):
    """Base class for an orbital, which can be evaluated at any point in space."""

    model_config = ConfigDict(arbitrary_types_allowed=True)  # Forbid extra fields not defined in the model
    type: str
    atom: str | None = None
    coords: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @cached_property
    def contraction_coeffs(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([])
    
    @cached_property
    def exponent_coeffs(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([])
    
    @cached_property
    def normalization_constants(self) -> np.ndarray:
        """Calculate the normalization constant for the orbital."""
        return np.array([])

    @abstractmethod
    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the orbital at a given point (x, y, z)."""
        pass


class STO3GOrbital(Orbital):
    """STO-3G for an orbital (1s, 2s, 2p, etc.), which is a linear combination of 3 GTOs."""
    # Exponent and coefficient for GTOs in the STO-3G expansion
    params1: tuple[float, float, int, int, int]  # (c1, alpha1, nx1, ny1, nz1)
    params2: tuple[float, float, int, int, int]  # (c2, alpha2, nx2, ny2, nz2)
    params3: tuple[float, float, int, int, int]  # (c3, alpha3, nx3, ny3, nz3)

    @cached_property
    def contraction_coeffs(self) -> np.ndarray:
        """Return the contraction coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[0], self.params2[0], self.params3[0]])

    @cached_property
    def exponent_coeffs(self) -> np.ndarray:
        """Return the exponent coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[1], self.params2[1], self.params3[1]])

    @cached_property
    def normalization_constants(self) -> np.ndarray:
        """Calculate the normalization constant for the STO-3G orbital."""
        return np.array([self.GTO1.N_cartesian, self.GTO2.N_cartesian, self.GTO3.N_cartesian])

    @cached_property
    def GTO1(self) -> GTO:
        return GTO(coords=self.coords, alpha=self.params1[1], nx=self.params1[2], ny=self.params1[3], nz=self.params1[4])
    
    @cached_property
    def GTO2(self) -> GTO:
        return GTO(coords=self.coords, alpha=self.params2[1], nx=self.params2[2], ny=self.params2[3], nz=self.params2[4])
    
    @cached_property
    def GTO3(self) -> GTO:
        return GTO(coords=self.coords, alpha=self.params3[1], nx=self.params3[2], ny=self.params3[3], nz=self.params3[4])

    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the STO-3G orbital at a given point (x, y, z)."""
        
        gto1 = self.GTO1(x, y, z) * self.params1[0]
        gto2 = self.GTO2(x, y, z) * self.params2[0]
        gto3 = self.GTO3(x, y, z) * self.params3[0]

        return gto1 + gto2 + gto3 # Use cartesian coordinates