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
from scipy.special import hyp1f1


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
        output = float(np.sum(matrix).item())
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
        output = float(np.sum(matrix).item())
        return 0.0 if math.isclose(output, 0.0) else output
    
    def VijR(self, orb1: STOGOrbital, orb2: STOGOrbital, R: tuple[float, float, float]) -> float:
        """Calculate the electron-nuclear attraction integral between two STO-nG orbitals for a single nucleus."""

        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u, cc1 in enumerate(orb1.cc):
            for v, cc2 in enumerate(orb2.cc):
                gto1 = orb1.gtos[u]
                gto2 = orb2.gtos[v]
                matrix[u, v] = gto1.N*gto2.N*cc1*cc2*self.compute_nuclear_attraction(gto1, gto2, R)
        output = float(np.sum(matrix).item())
        return 0.0 if math.isclose(output, 0.0) else output


    def Vijkl(self, orb1: STOGOrbital, orb2: STOGOrbital, orb3: STOGOrbital, orb4: STOGOrbital) -> float:
        """Calculate the Coulomb repulsion integral between two STO-nG orbitals."""
        
        matrix = np.zeros((len(orb1.cc), len(orb2.cc), len(orb3.cc), len(orb4.cc)))
        for m, cc1 in enumerate(orb1.cc):
            for n, cc2 in enumerate(orb2.cc):
                for u, cc3 in enumerate(orb3.cc):
                    for v, cc4 in enumerate(orb4.cc):
                        gto1 = orb1.gtos[m]
                        gto2 = orb2.gtos[n]
                        gto3 = orb3.gtos[u]
                        gto4 = orb4.gtos[v]

                        norms = gto1.N * gto2.N * gto3.N * gto4.N
                        coefs = cc1 * cc2 * cc3 * cc4
                        matrix[m, n, u, v] = norms * coefs * self.compute_repulsion(gto1, gto2, gto3, gto4)
        output = float(np.sum(matrix).item()) + 1 - 1
        return 0.0 if math.isclose(output, 0.0) else output
                        
    def compute_repulsion(self, gto1: STOGPrimitive, gto2: STOGPrimitive, gto3: STOGPrimitive, gto4: STOGPrimitive) -> float:
        """Calculate the electron-electron repulsion integral between two STO-nG orbitals. https://joshuagoings.com/2017/04/28/integrals/"""
        p = gto1.alpha + gto2.alpha
        q = gto3.alpha + gto4.alpha
        alpha = p * q / (p + q)
        P = [self.compute_Pi(gto1, gto2, component) for component in ['x', 'y', 'z']]
        Q = [self.compute_Pi(gto3, gto4, component) for component in ['x', 'y', 'z']]

        val = 0.0
        for t in range(gto1.nx+gto2.nx+1):
            for u in range(gto1.ny+gto2.ny+1):
                for v in range(gto1.nz+gto2.nz+1):
                    for tau in range(gto3.nx+gto4.nx+1):
                        for nu in range(gto3.ny+gto4.ny+1):
                            for phi in range(gto3.nz+gto4.nz+1):
                                term1 = self.compute_E(gto1.nx, gto2.nx, t, gto1.coords[0], gto2.coords[0], gto1.alpha, gto2.alpha)
                                term2 = self.compute_E(gto1.ny, gto2.ny, u, gto1.coords[1], gto2.coords[1], gto1.alpha, gto2.alpha)
                                term3 = self.compute_E(gto1.nz, gto2.nz, v, gto1.coords[2], gto2.coords[2], gto1.alpha, gto2.alpha)
                                term4 = self.compute_E(gto3.nx, gto4.nx, tau, gto3.coords[0], gto4.coords[0], gto3.alpha, gto4.alpha)
                                term5 = self.compute_E(gto3.ny, gto4.ny, nu, gto3.coords[1], gto4.coords[1], gto3.alpha, gto4.alpha)
                                term6 = self.compute_E(gto3.nz, gto4.nz, phi, gto3.coords[2], gto4.coords[2], gto3.alpha, gto4.alpha)
                                term7 = np.power(-1, tau+nu+phi)
                                term8 = self.compute_R(t+tau, u+nu, v+phi, 0, alpha, P, Q)
                                val += term1*term2*term3*term4*term5*term6*term7*term8
        prefactor = 2 * math.pi**(5/2) / (p * q * math.sqrt(p + q))
        return prefactor * val

    def compute_nuclear_attraction(self, gto1: STOGPrimitive, gto2: STOGPrimitive, R: tuple[float, float, float]) -> float:
        """Calculate the nuclear attraction integral between two STO-nG orbitals for a single nucleus."""
        p = gto1.alpha + gto2.alpha
        P = [self.compute_Pi(gto1, gto2, component) for component in ['x', 'y', 'z']] # Gaussian composite center
        val = 0.0
        for t in range(gto1.nx+gto2.nx+1):
            for u in range(gto1.ny+gto2.ny+1):
                for v in range(gto1.nz+gto2.nz+1):
                    val += self.compute_E(gto1.nx, gto2.nx, t, gto1.coords[0], gto2.coords[0], gto1.alpha, gto2.alpha) * \
                        self.compute_E(gto1.ny, gto2.ny, u,gto1.coords[1], gto2.coords[1], gto1.alpha, gto2.alpha) * \
                        self.compute_E(gto1.nz, gto2.nz, v,gto1.coords[2], gto2.coords[2], gto1.alpha, gto2.alpha) * \
                        self.compute_R(t, u, v, 0, p, P, R)
        return 2*np.pi/p  * val

    def compute_tx(self, A: float, B: float, alpha: float, beta: float, ai: int, bi: int) -> float:
        """Calculate the integal between 2 Gaussion primitives in a given direction for the kinetic energy integral.
        https://content.wolfram.com/sites/19/2013/01/Ho_Kinetic.pdf"""

        P = (alpha*A + beta*B) / (alpha+beta)
        si = partial(self.compute_sx, A=A, B=B, alpha=alpha, beta=beta)
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

    ####### Common subroutines for Sij, Tij, and Vij #######

    def compute_E(self, ai: int, bi: int, t: int, A: float, B: float, alpha: float, beta: float) -> float:
        """ Recursive definition of Hermite Gaussian coefficients.
            a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
            b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
            i,j: orbital angular momentum number on Gaussian 'a' and 'b'
            t: number nodes in Hermite (depends on type of integral, e.g. always zero for overlap integrals)
            A: origin of Gaussian 'a'
            B: origin of Gaussian 'b'
        """
        p = alpha + beta
        q = alpha*beta/p
        Qx = A - B
        if (t < 0) or (t > (ai + bi)):
            return 0.0
        if ai == bi == t == 0:
            return np.exp(-q*Qx**2)
        if bi == 0:
            # decrement index i
            return (1/(2*p))*self.compute_E(ai-1,bi,t-1,A,B,alpha,beta) - \
                (q*Qx/alpha)*self.compute_E(ai-1,bi,t,A,B,alpha,beta)    + \
                (t+1)*self.compute_E(ai-1,bi,t+1,A,B,alpha,beta)
        # decrement index j
        return (1/(2*p))*self.compute_E(ai,bi-1,t-1,A,B,alpha,beta) + \
            (q*Qx/beta)*self.compute_E(ai,bi-1,t,A,B,alpha,beta)    + \
            (t+1)*self.compute_E(ai,bi-1,t+1,A,B,alpha,beta)

    def boys(self, n: int, T: float) -> float:
        return hyp1f1(n+0.5,n+1.5,-T)/(2.0*n+1.0) 

    def compute_R(self, t: int,u: int, v:int, n:int, p:float, P: tuple[float, float, float], C: tuple[float, float, float]) -> float:
        ''' Returns the Coulomb auxiliary Hermite integrals 
            Arguments:
                t,u,v: order of Coulomb Hermite derivative in x,y,z
                n: order of Boys function 
                PC: Gaussian composite center P.
                C: Nuclear center C.
        '''
        RPC = math.sqrt((P[0]-C[0])**2 + (P[1]-C[1])**2 + (P[2]-C[2])**2)
        T = p*RPC*RPC
        val = 0.0
        if t == u == v == 0:
            val += np.power(-2*p,n)*self.boys(n,T)
        elif t == u == 0:
            if v > 1:
                val += (v-1)*self.compute_R(t,u,v-2,n+1,p,P,C)
            val += (P[2]-C[2])*self.compute_R(t,u,v-1,n+1,p,P,C)
        elif t == 0:
            if u > 1:
                val += (u-1)*self.compute_R(t,u-2,v,n+1,p,P,C)
            val += (P[1]-C[1])*self.compute_R(t,u-1,v,n+1,p,P,C)
        else:
            if t > 1:
                val += (t-1)*self.compute_R(t-2,u,v,n+1,p,P,C)
            val += (P[0]-C[0])*self.compute_R(t-1,u,v,n+1,p,P,C)
        return val

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
        
    def compute_sx(self, A: float, B: float, alpha: float, beta: float, ai: int, bi: int) -> float:
        """Calculate the integral of the product of two GTOs in a given component."""
        P = (alpha*A + beta*B) / (alpha+beta)

        if ai < 0 or bi < 0:
            return 0
        if (ai, bi) == (0, 0):
            return 1
        if (ai, bi) == (1, 0):
            return -(A-P)
        
        func = partial(self.compute_sx, A=A, B=B, alpha=alpha, beta=beta)
        if ai > 1 and bi == 0:
            return -(A-P)*func(ai=ai-1, bi=0) + ((ai-1)/(2*alpha + 2*beta)) * func(ai=ai-2, bi=0)
        return func(ai=ai+1, bi=bi-1) + (A-B)*func(ai=ai, bi=bi-1)
    
    def compute_integral_by_func(self, gto1: STOGPrimitive, gto2: STOGPrimitive, func: Callable, component: Literal['x', 'y', 'z'] = 'x', R: tuple[float, float, float] = None, **kwargs) -> float:
        """Compute the integral between two GTO primitives in a given component using the specified function."""
        A = getattr(gto1, f"coords")[['x', 'y', 'z'].index(component)]
        B = getattr(gto2, f"coords")[['x', 'y', 'z'].index(component)]
        ai = getattr(gto1, f"n{component}")
        bi = getattr(gto2, f"n{component}")
        if R is not None:
            Ri = R[['x', 'y', 'z'].index(component)]
            return func(A, B, gto1.alpha, gto2.alpha, ai, bi, R=Ri, **kwargs)
        return func(A, B, gto1.alpha, gto2.alpha, ai, bi, **kwargs)
    
    def get_factors(self, orb1: STOGOrbital, orb2: STOGOrbital, u: int, v: int) -> tuple[float, float, float]:
        """Calculate the common factors for the integrals between two GTO primitives."""
        coeff = orb1.cc[u] * orb1.N[u] * orb2.cc[v] * orb2.N[v]
        E_AB = self.compute_E_AB(orb1.gtos[u], orb2.gtos[v])
        prefactor = (math.pi / (orb1.gtos[u].alpha + orb2.gtos[v].alpha)) ** (3 / 2)
        return coeff, E_AB, prefactor
