import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field
from functools import cached_property
from typing import Callable, Literal
# from source.hfp import GTO_LIB
import numpy as np
from scipy.special import factorial2
from functools import partial
from hf.chebyshev import IntChebyshev



class STOGPrimitive(BaseModel):
    """Gaussian Type Orbital (GTO) basis function."""
    coords: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)

    # Exponent and angular momentum quantum numbers for the GTO (s, p, d, etc.)
    alpha: float
    nx: int = 0
    ny: int = 0
    nz: int = 0

    @property
    def l(self) -> int:
        """Calculate the total angular momentum quantum number."""
        return self.nx + self.ny + self.nz

    @cached_property
    def N(self) -> float:
        """Normalization constant for the cartisian GTO. https://en.wikipedia.org/wiki/Gaussian_orbital#Cartesian_coordinates"""
        prefactor = (2 * self.alpha / math.pi) ** (3 / 4)
        numerator = (8 * self.alpha) ** (self.nx + self.ny + self.nz) * math.factorial(self.nx) * math.factorial(self.ny) * math.factorial(self.nz)
        denominator = math.factorial(2*self.nx) * math.factorial(2*self.ny) * math.factorial(2*self.nz)
        return prefactor * math.sqrt(numerator/ denominator)
    

    def Y(self, x: float, y: float , z: float) -> float:
        """Calculate the angular part of the GTO in cartesian coordinates."""
        return (x - self.coords[0]) ** self.nx * (y - self.coords[1]) ** self.ny * (z - self.coords[2]) ** self.nz

    def r2(self, x: float, y: float, z: float) -> float:
        """Calculate the squared distance from the center of the GTO."""
        dist_x = x - self.coords[0]
        dist_y = y - self.coords[1]
        dist_z = z - self.coords[2]
        return dist_x ** 2 + dist_y ** 2 + dist_z ** 2

    def r(self, x: float, y: float, z: float) -> float:
        """Calculate the distance from the center of the GTO."""
        return math.sqrt(self.r2(x, y, z))

    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the GTO at a given point (x, y, z)."""
        R = math.exp(-self.alpha * self.r2(x, y, z))
        return self.N * R * self.Y(x, y, z) # Use cartesian coordinates

class STOGOrbital(BaseModel):
    """STO-nG for an atomic orbital (1s, 2s, 2p, etc.), which is a linear combination of n GTOs."""

    model_config = ConfigDict(arbitrary_types_allowed=True)  # Forbid extra fields not defined in the model
    
    type: str
    atom: str | None = None
    coords: tuple[float, float, float] = (0.0, 0.0, 0.0)

    cc: list[float] = Field(..., description="Contraction coefficients for the GTOs in the STO-nG expansion")
    alpha: list[float] = Field(..., description="Exponent coefficients for the GTOs in the STO-nG expansion")
    nx: list[int] = Field(..., description="Angular momentum quantum numbers in x direction for the GTOs in the STO-nG expansion")
    ny: list[int] = Field(..., description="Angular momentum quantum numbers in y direction for the GTOs in the STO-nG expansion")
    nz: list[int] = Field(..., description="Angular momentum quantum numbers in z direction for the GTOs in the STO-nG expansion")

    @property
    def l(self) -> int:
        """Calculate the total angular momentum quantum number."""
        return self.nx[0] + self.ny[0] + self.nz[0]

    @cached_property
    def gtos(self) -> list[STOGPrimitive]:
        """Calculate the GTOs that make up the STO-nG orbital."""
        return [STOGPrimitive(coords=self.coords, alpha=self.alpha[i], nx=self.nx[i], ny=self.ny[i], nz=self.nz[i]) for i in range(len(self.cc))]
    
    @cached_property
    def N(self) -> float:
        """The normalization constant for the STO-nG orbital."""
        return [gto.N for gto in self.gtos]
    
    def __call__(self, x: float, y: float, z: float) -> float:
        """Evaluate the STO-nG orbital at a given point (x, y, z)."""
        return sum(self.cc[i] * self.gtos[i](x, y, z) for i in range(len(self.cc))) # Use cartesian coordinates

class STOGIntegrator:
    """
    Integrator class for calculating integrals between STO-nG orbitals.
    The idea is for each S(i,j) where i is the ith atomic orbital, we compute
    a mxn matrix where each element is an integral between Gausian primitive
    m of atomic orbital i and Gaussian primitive n of atomic orbital j.
    We than sum over all the elements in the matrix to get S(i,j).
    We do this for all pairs of atomic orbitals to get the full overlap matrix S.
    """

    def __init__(self):
        self.cheby_shev_integrator = IntChebyshev()

    def Sij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the overlap integral between two STO-nG orbitals."""
    
        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u in range(n):
            for v in range(m):
                coeff, E_AB, prefactor = self.get_factors(orb1, orb2, u, v)
                s_x = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'x')
                s_y = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'y')
                s_z = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'z')
                matrix[u, v] = coeff * E_AB * prefactor * s_x * s_y * s_z
        output = float(np.sum(matrix).item()) + 1 - 1
        return 0.0 if math.isclose(output, 0.0) else output
    
    def Tij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the kinetic energy integral between two STO-nG orbitals."""
    
        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u in range(n):
            for v in range(m):
                coeff, E_AB, prefactor = self.get_factors(orb1, orb2, u, v)
                t_x = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_tx, 'x')
                t_y = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_tx, 'y')
                t_z = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_tx, 'z')
                s_x = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'x')
                s_y = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'y')
                s_z = self.compute_integral_by_func(orb1.gtos[u], orb2.gtos[v], self.compute_sx, 'z')
                matrix[u, v] = coeff * E_AB * prefactor * (t_x*s_y*s_z + t_y*s_x*s_z + t_z*s_x*s_y) 
        output = float(np.sum(matrix).item()) + 1 - 1
        return 0.0 if math.isclose(output, 0.0) else output
    
    def VijR(self, orb1: STOGOrbital, orb2: STOGOrbital, R: tuple[float, float, float]) -> float:
        """Calculate the electron-nuclear attraction integral between two STO-nG orbitals for a single nucleus."""
        
        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u in range(n):
            for v in range(m):
                coeff, E_AB, _ = self.get_factors(orb1, orb2, u, v)
                gto1 = orb1.gtos[u]
                gto2 = orb2.gtos[v]

                def func(t: float) -> float:
                    term1 = (gto1.alpha + gto2.alpha)*t**2
                    term2 = (gto1.alpha * np.array(gto1.coords) + gto2.alpha * np.array(gto2.coords)) / (gto1.alpha + gto2.alpha) - np.array(R)
                    term3 = np.dot(term2, term2)
                    v_x = self.compute_integral_by_func(gto1, gto2, self.compute_nx, 'x', R=R, t=t)
                    v_y = self.compute_integral_by_func(gto1, gto2, self.compute_nx, 'y', R=R, t=t)
                    v_z = self.compute_integral_by_func(gto1, gto2, self.compute_nx, 'z', R=R, t=t)
                    return math.exp(-term1*term3) * v_x * v_y * v_z
                prefactor = 2 * math.pi / (gto1.alpha + gto2.alpha)
                matrix[u, v] = coeff * E_AB * prefactor * self.cheby_shev_integrator.integrate(eps=1e-10, f=lambda x: 0.5 * func(0.5*(x+1)), m=50000)
        output = float(np.sum(matrix).item()) + 1 - 1
        return 0.0 if math.isclose(output, 0.0) else output

    def Jij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the Coulomb integral between two STO-nG orbitals."""
        pass

    def Kij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the exchange integral between two STO-nG orbitals."""
        pass

    def compute_tx(self, A: float, B: float, P: float, alpha: float, beta: float, ai: int, bi: int) -> float:
        """Calculate the integal between 2 Gaussion primitives in a given direction for the kinetic energy integral.
        https://content.wolfram.com/sites/19/2013/01/Ho_Kinetic.pdf"""

        si = partial(self.compute_sx, A=A, B=B, P=P, alpha=alpha, beta=beta)
        if ai < 0 or bi < 0:
            return 0
        if alpha == 0 and beta == 0:
            return 2*alpha*beta*si(ai=1,bi=1)
        if bi == 0:
            return -ai*beta*si(ai=ai-1,bi=1) + 2*alpha*beta*si(ai=ai+1,bi=1)
        if ai == 0:
            return -alpha*bi*si(ai=1,bi=bi-1) + 2*alpha*beta*si(ai=1,bi=bi+1)
        term1 = ai*bi*si(ai=ai-1,bi=bi-1)
        term2 = 2*ai*beta*si(ai=ai-1,bi=bi+1)
        term3 = 2*alpha*bi*si(ai=ai+1,bi=bi-1)
        term4 = 4*alpha*beta*si(ai=ai+1,bi=bi+1)
        return 0.5*(term1-term2-term3+term4)

    def compute_nx(self, A: float, B: float, P: float, alpha: float, beta: float, ai: int, bi: int, t: float, R: float) -> float:
        """https://content.wolfram.com/sites/19/2014/12/Ho_Nuclear.pdf"""
        func = partial(self.compute_nx, A=A, B=B, P=P, alpha=alpha, beta=beta, t=t, R=R)
        if ai < 0 or bi < 0:
            return 0
        if (ai, bi) == (0, 0):
            return 1
        if (ai, bi) == (1, 0):
            result = -(A-P) - t**2*(P-R)
            return result
        if bi == 0:
            term1 = -(A-P)
            term2 = -t**2 * (P-R) * func(ai=ai-1, bi=0)
            term3 = (ai-1)/(2*alpha + 2*beta) * (1-t**2) * func(ai=ai-2, bi=0)
            return term1 + term2 + term3
        return func(ai=ai+1, bi=bi-1) + (A-B)*func(ai=ai, bi=bi-1)

    ####### Common subroutines for Sij, Tij, and Vij #######
    def compute_E_AB(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the exponential prefactor E_AB for the overlap integral."""
        prefactor = - gto1.alpha * gto2.alpha / (gto1.alpha + gto2.alpha)
        r_squared = self.compute_r2(gto1, gto2)
        return math.exp(prefactor * r_squared)
    
    def compute_r2(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the squared distance between the centers of two GTO primitives."""
        dist_x = gto1.coords[0] - gto2.coords[0]
        dist_y = gto1.coords[1] - gto2.coords[1]
        dist_z = gto1.coords[2] - gto2.coords[2]
        return dist_x ** 2 + dist_y ** 2 + dist_z ** 2
    
    def compute_Pi(self, gto1: STOGPrimitive, gto2: STOGPrimitive, component: Literal['x', 'y', 'z'] = 'x') -> float:
        """Calculate the weighted center and combined exponent for the integral in a given component."""
        alpha_sum = gto1.alpha + gto2.alpha
        if component == 'x':
            return (gto1.alpha * gto1.coords[0] + gto2.alpha * gto2.coords[0]) / alpha_sum
        elif component == 'y':
            return (gto1.alpha * gto1.coords[1] + gto2.alpha * gto2.coords[1]) / alpha_sum
        else:  # component == 'z'
            return (gto1.alpha * gto1.coords[2] + gto2.alpha * gto2.coords[2]) / alpha_sum
        
    def compute_sx(self, A: float, B: float, P: float, alpha: float, beta: float, ai: int, bi: int) -> float:
        """Calculate the integral of the product of two GTOs in a given component."""

        if ai < 0 or bi < 0:
            return 0
        if (ai, bi) == (0, 0):
            return 1
        if (ai, bi) == (1, 0):
            return -(A-P)
        
        func = partial(self.compute_sx, A=A, B=B, P=P, alpha=alpha, beta=beta)
        if ai > 1 and bi == 0:
            return -(A-P)*func(ai=ai-1, bi=0) + ((ai-1)/(2*alpha + 2*beta)) * func(ai=ai-2, bi=0)
        return func(ai=ai+1, bi=bi-1) + (A-B)*func(ai=ai, bi=bi-1)
    
    def compute_integral_by_func(self, gto1: STOGPrimitive, gto2: STOGPrimitive, func: Callable, component: Literal['x', 'y', 'z'] = 'x', R: tuple[float, float, float] = None, **kwargs) -> float:
        """Compute the integral between two GTO primitives in a given component using the specified function."""
        A = getattr(gto1, f"coords")[['x', 'y', 'z'].index(component)]
        B = getattr(gto2, f"coords")[['x', 'y', 'z'].index(component)]
        P = self.compute_Pi(gto1, gto2, component)
        ai = getattr(gto1, f"n{component}")
        bi = getattr(gto2, f"n{component}")
        if R is not None:
            Ri = R[['x', 'y', 'z'].index(component)]
            return func(A, B, P, gto1.alpha, gto2.alpha, ai, bi, R=Ri, **kwargs)
        return func(A, B, P, gto1.alpha, gto2.alpha, ai, bi, **kwargs)
    
    def get_factors(self, orb1: STOGOrbital, orb2: STOGOrbital, u: int, v: int) -> tuple[float, float, float]:
        """Calculate the common factors for the integrals between two GTO primitives."""
        coeff = orb1.cc[u] * orb1.N[u] * orb2.cc[v] * orb2.N[v]
        E_AB = self.compute_E_AB(orb1.gtos[u], orb2.gtos[v])
        prefactor = (math.pi / (orb1.gtos[u].alpha + orb2.gtos[v].alpha)) ** (3 / 2)
        return coeff, E_AB, prefactor