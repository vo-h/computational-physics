import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field, computed_field
from functools import cached_property
import numpy as np
from abc import ABC, abstractmethod

class GTO(BaseModel):
    """Gaussian Type Orbital (GTO) basis function."""
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)

    # Exponent and angular momentum quantum numbers for the GTO (s, p, d, etc.)
    alpha: float
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

    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the GTO at a given point (x, y, z)."""
        r_squared = (x - self.coordinates[0]) ** 2 + (y - self.coordinates[1]) ** 2 + (z - self.coordinates[2]) ** 2
        return self.N * (x - self.coordinates[0]) ** self.nx * (y - self.coordinates[1]) ** self.ny * (z - self.coordinates[2]) ** self.nz * math.exp(-self.alpha * r_squared)


class Orbital(BaseModel, ABC):
    """Base class for an orbital, which can be evaluated at any point in space."""

    model_config = ConfigDict(arbitrary_types_allowed=True)  # Forbid extra fields not defined in the model
    type: str
    atom: str | None = None
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @cached_property
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([])
    
    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
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
    def prime_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[0], self.params2[0], self.params3[0]])

    @cached_property
    def orbital_coefficients(self) -> np.ndarray:
        """Return the coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.array([self.params1[1], self.params2[1], self.params3[1]])

    @cached_property
    def normalization_constants(self) -> np.ndarray:
        """Calculate the normalization constant for the STO-3G orbital."""
        return np.array([self.GTO1.N, self.GTO2.N, self.GTO3.N])

    @cached_property
    def GTO1(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params1[1], nx=self.params1[2], ny=self.params1[3], nz=self.params1[4])
    
    @cached_property
    def GTO2(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params2[1], nx=self.params2[2], ny=self.params2[3], nz=self.params2[4])
    
    @cached_property
    def GTO3(self) -> GTO:
        return GTO(coordinates=self.coordinates, alpha=self.params3[1], nx=self.params3[2], ny=self.params3[3], nz=self.params3[4])

    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the STO-3G orbital at a given point (x, y, z)."""
        gto1 = self.GTO1(x, y, z) * self.params1[0]
        gto2 = self.GTO2(x, y, z) * self.params2[0]
        gto3 = self.GTO3(x, y, z) * self.params3[0]
        return gto1 + gto2 + gto3