import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field
from functools import cached_property
from typing import Callable, Literal
# from source.hfp import GTO_LIB
import numpy as np
from scipy.special import factorial2
from functools import partial
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

    def Sij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the overlap integral between two STO-nG orbitals."""

        def compute_Sab(gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
            """Calculate the overlap integral S_ab between two GTO primitives. (eq.101)"""
            prefactor = (math.pi / (gto1.alpha + gto2.alpha)) ** (3 / 2)
            s_x = self.compute_E(gto1.coords[0], gto2.coords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, 0)
            s_y = self.compute_E(gto1.coords[1], gto2.coords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, 0)
            s_z = self.compute_E(gto1.coords[2], gto2.coords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, 0)
            return prefactor * s_x * s_y * s_z

        sij = 0.0
        for u, cc1 in enumerate(orb1.cc):
            for v, cc2 in enumerate(orb2.cc):
                sij += cc1 * cc2 * orb1.gtos[u].N * orb2.gtos[v].N * compute_Sab(orb1.gtos[u], orb2.gtos[v])
        return 0.0 if math.isclose(sij, 0.0) else sij
    
    def Tij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the kinetic energy integral between two STO-nG orbitals."""

        def compute_D(A: float, B: float, alpha: float, beta: float, ai: int, bi: int, t: int) -> float:
            """ Returns the Dawson function D_n(T) """
            p = alpha + beta
            if t == 0: # eq. 113 
                return self.compute_E(A, B, alpha, beta, ai, bi, 0) * np.sqrt(np.pi/(p))
            term1 = bi * compute_D(A, B, alpha, beta, ai, bi-1, t-1)
            term2 = 2*beta*compute_D(A, B, alpha, beta, ai, bi+1, t-1)
            return term1 - term2 # eq. 114

        def compute_Tab(orb1: STOGOrbital, orb2: STOGOrbital) -> float:
            """Calculate the kinetic energy integral T_ab between two GTO primitives. (eq. 121)"""
            term1 = compute_D(orb1.coords[0], orb2.coords[0], orb1.alpha, orb2.alpha, orb1.nx, orb2.nx, 2)
            term2 = compute_D(orb1.coords[1], orb2.coords[1], orb1.alpha, orb2.alpha, orb1.ny, orb2.ny, 0)
            term3 = compute_D(orb1.coords[2], orb2.coords[2], orb1.alpha, orb2.alpha, orb1.nz, orb2.nz, 0)
            term4 = compute_D(orb1.coords[0], orb2.coords[0], orb1.alpha, orb2.alpha, orb1.nx, orb2.nx, 0)
            term5 = compute_D(orb1.coords[1], orb2.coords[1], orb1.alpha, orb2.alpha, orb1.ny, orb2.ny, 2)
            term6 = compute_D(orb1.coords[2], orb2.coords[2], orb1.alpha, orb2.alpha, orb1.nz, orb2.nz, 0)
            term7 = compute_D(orb1.coords[0], orb2.coords[0], orb1.alpha, orb2.alpha, orb1.nx, orb2.nx, 0)
            term8 = compute_D(orb1.coords[1], orb2.coords[1], orb1.alpha, orb2.alpha, orb1.ny, orb2.ny, 0)
            term9 = compute_D(orb1.coords[2], orb2.coords[2], orb1.alpha, orb2.alpha, orb1.nz, orb2.nz, 2)
            return -0.5 * (term1*term2*term3 + term4*term5*term6 + term7*term8*term9)

        tij = 0.0
        for u, cc1 in enumerate(orb1.cc):
            for v, cc2 in enumerate(orb2.cc):
                tij += cc1 * cc2 * orb1.gtos[u].N * orb2.gtos[v].N * compute_Tab(orb1.gtos[u], orb2.gtos[v])
        return 0.0 if math.isclose(tij, 0.0) else tij
    
    def VijR(self, orb1: STOGOrbital, orb2: STOGOrbital, R: tuple[float, float, float]) -> float:
        """Calculate the electron-nuclear attraction integral between two STO-nG orbitals for a single nucleus."""

        def compute_Vab(gto1: STOGPrimitive, gto2: STOGPrimitive, R: tuple[float, float, float]) -> float:
            """Calculate the nuclear attraction integral between two STO-nG orbitals for a single nucleus. (eq. 199-204)"""
            p = gto1.alpha + gto2.alpha
            P = [self.compute_Pi(gto1, gto2, component) for component in ['x', 'y', 'z']] # Gaussian composite center
            val = 0.0
            for t in range(gto1.nx+gto2.nx+1):
                for u in range(gto1.ny+gto2.ny+1):
                    for v in range(gto1.nz+gto2.nz+1):
                        val += self.compute_E(gto1.coords[0], gto2.coords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t) * \
                            self.compute_E(gto1.coords[1], gto2.coords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u) * \
                            self.compute_E(gto1.coords[2], gto2.coords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v) * \
                            self.compute_R(t, u, v, 0, p, P, R)
            return 2*np.pi/p  * val

        vijr = 0.0
        for u, cc1 in enumerate(orb1.cc):
            for v, cc2 in enumerate(orb2.cc):
                vijr += orb1.gtos[u].N*orb2.gtos[v].N*cc1*cc2*compute_Vab(orb1.gtos[u], orb2.gtos[v], R)
        return 0.0 if math.isclose(vijr, 0.0) else vijr

    def Vijkl(self, orb1: STOGOrbital, orb2: STOGOrbital, orb3: STOGOrbital, orb4: STOGOrbital) -> float:
        """Calculate the Coulomb repulsion integral between two STO-nG orbitals."""

        def compute_Vabcd(gto1: STOGPrimitive, gto2: STOGPrimitive, gto3: STOGPrimitive, gto4: STOGPrimitive) -> float:
            """Calculate the electron-electron repulsion integral between two STO-nG orbitals."""
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
                                    term1 = self.compute_E(gto1.coords[0], gto2.coords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t)
                                    term2 = self.compute_E(gto1.coords[1], gto2.coords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u)
                                    term3 = self.compute_E(gto1.coords[2], gto2.coords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v)
                                    term4 = self.compute_E(gto3.coords[0], gto4.coords[0], gto3.alpha, gto4.alpha, gto3.nx, gto4.nx, tau)
                                    term5 = self.compute_E(gto3.coords[1], gto4.coords[1], gto3.alpha, gto4.alpha, gto3.ny, gto4.ny, nu)
                                    term6 = self.compute_E(gto3.coords[2], gto4.coords[2], gto3.alpha, gto4.alpha, gto3.nz, gto4.nz, phi)
                                    term7 = np.power(-1, tau+nu+phi)
                                    term8 = self.compute_R(t+tau, u+nu, v+phi, 0, alpha, P, Q)
                                    val += term1*term2*term3*term4*term5*term6*term7*term8
            prefactor = 2 * math.pi**(5/2) / (p * q * math.sqrt(p + q))
            return prefactor * val
        
        vijkl = 0.0
        for m, cc1 in enumerate(orb1.cc):
            for n, cc2 in enumerate(orb2.cc):
                for u, cc3 in enumerate(orb3.cc):
                    for v, cc4 in enumerate(orb4.cc):
                        norms = orb1.gtos[m].N * orb2.gtos[n].N * orb3.gtos[u].N * orb4.gtos[v].N
                        coefs = cc1 * cc2 * cc3 * cc4
                        vijkl += norms * coefs * compute_Vabcd(orb1.gtos[m], orb2.gtos[n], orb3.gtos[u], orb4.gtos[v])
        return 0.0 if math.isclose(vijkl, 0.0) else vijkl


    ####### Common subroutines for Sij, Tij, and Vij #######

    def compute_E(self, A: float, B: float, alpha: float, beta: float, ai: int, bi: int, t: int) -> float:
        """ Recursive definition of Hermite Gaussian coefficients using eq. 73-75.
            A, B: origin of Gaussian 'a' and 'b'
            alpha, beta: exponent of Gaussian 'a' and 'b'
            ai, bi: orbital angular momentum number on Gaussian 'a' and 'b'. i,j in the paper; nx, ny, nz in the code.
            t: number nodes in Hermite (depends on type of integral, e.g. always zero for overlap integrals)
        """

        func = partial(self.compute_E, A=A, B=B, alpha=alpha, beta=beta)
        p = alpha + beta
        q = alpha*beta/p
        Qx = A - B
        if (t < 0) or (t > (ai + bi)):
            return 0.0
        if ai == bi == t == 0: # eq: 73
            return np.exp(-q*Qx**2)
        if bi == 0: # eq: 74
            term1 = (1/(2*p))*func(ai=ai-1,bi=bi,t=t-1)
            term2 = (q*Qx/alpha)*func(ai=ai-1,bi=bi,t=t)
            term3 = (t+1)*func(ai=ai-1,bi=bi,t=t+1)
            return term1 - term2 + term3
        # eq: 75
        term1 = (1/(2*p))*func(ai=ai,bi=bi-1,t=t-1)
        term2 = (q*Qx/beta)*func(ai=ai,bi=bi-1,t=t)
        term3 = (t+1)*func(ai=ai,bi=bi-1,t=t+1)
        return term1 + term2 + term3

    def compute_R(self, t: int,u: int, v:int, n:int, p:float, P: tuple[float, float, float], C: tuple[float, float, float]) -> float:
        """ Returns the Coulomb auxiliary Hermite integrals (eq. 190-192)
            Arguments:
                t,u,v: order of Coulomb Hermite derivative in x,y,z
                n: order of Boys function 
                C: Nuclear center C.
        """

        def boys(n: int, T: float) -> float:
            """Boys function using confluent hypergeometric function."""
            return hyp1f1(n+0.5,n+1.5,-T)/(2.0*n+1.0) 

        RPC = math.sqrt((P[0]-C[0])**2 + (P[1]-C[1])**2 + (P[2]-C[2])**2)
        val = 0.0
        if t < 0 or u < 0 or v < 0:
            return 0.0
        if t == u == v == 0: # eq. 182
            val += (-2*p)**(n)*boys(n,p*RPC*RPC)
        elif t == u == 0:
            val +=(v-1)*self.compute_R(t,u,v-2,n+1,p,P,C)
            val += (P[2]-C[2])*self.compute_R(t,u,v-1,n+1,p,P,C)
        elif t == 0:
            val += (u-1)*self.compute_R(t,u-2,v,n+1,p,P,C)
            val += (P[1]-C[1])*self.compute_R(t,u-1,v,n+1,p,P,C)
        else:
            val += (t-1)*self.compute_R(t-2,u,v,n+1,p,P,C)
            val += (P[0]-C[0])*self.compute_R(t-1,u,v,n+1,p,P,C)
        return val
    
    def compute_Pi(self, gto1: STOGPrimitive, gto2: STOGPrimitive, component: Literal['x', 'y', 'z'] = 'x') -> float:
        """Calculate the weighted center and combined exponent for the integral in a given component."""
        alpha_sum = gto1.alpha + gto2.alpha
        if component == 'x':
            return (gto1.alpha * gto1.coords[0] + gto2.alpha * gto2.coords[0]) / alpha_sum
        elif component == 'y':
            return (gto1.alpha * gto1.coords[1] + gto2.alpha * gto2.coords[1]) / alpha_sum
        else:  # component == 'z'
            return (gto1.alpha * gto1.coords[2] + gto2.alpha * gto2.coords[2]) / alpha_sum


    

