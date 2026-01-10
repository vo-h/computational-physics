import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field
from functools import cached_property
from typing import Callable, Literal
from source.hf import GTO_LIB
import numpy as np
from scipy.special import factorial2


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
    
    @cached_property
    def Y(self) -> Callable[[float, float, float], float]:
        """Calculate the angular part of the GTO in cartesian coordinates."""
        return lambda x, y, z: (x - self.coords[0]) ** self.nx * (y - self.coords[1]) ** self.ny * (z - self.coords[2]) ** self.nz

    def r2(self, x: float, y: float, z: float) -> float:
        """Calculate the squared distance from the center of the GTO."""
        dist_x = x - self.coords[0]
        dist_y = y - self.coords[1]
        dist_z = z - self.coords[2]
        return dist_x ** 2 + dist_y ** 2 + dist_z ** 2

    def r(self, x: float, y: float, z: float) -> float:
        """Calculate the distance from the center of the GTO."""
        return math.sqrt(self.r2(x, y, z))

    def __call__(self, x: float, y: float, z: float, resolver: Literal['C', 'Python'] = 'Python') -> float:
        """Evaluate the GTO at a given point (x, y, z)."""

        if resolver == 'Python':
            R = math.exp(-self.alpha * self.r2(x, y, z))
            return self.N * R * self.Y(x, y, z) # Use cartesian coordinates

        return GTO_LIB.gto(
            np.array(self.coords),
            np.array([x, y, z]),
            np.array([self.nx, self.ny, self.nz], dtype=np.int32),
            self.alpha
        )


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
    
    def __call__(self, x: float, y: float, z: float, resolver: Literal["C", "Python"] = "Python") -> float:
        """Evaluate the STO-nG orbital at a given point (x, y, z)."""
        if resolver == "Python":
            return sum(self.cc[i] * self.gtos[i](x, y, z) for i in range(len(self.cc))) # Use cartesian coordinates
        return GTO_LIB.sto_ng(
            np.array(self.coords),
            np.array([x, y, z]),
            np.array(self.cc),
            np.array(self.alpha),
            np.array(self.nx, dtype=np.int32),
            np.array(self.ny, dtype=np.int32),
            np.array(self.nz, dtype=np.int32),
            len(self.cc)
        )

class SijIntegrator:
    """Compute Sij for the overlap matrix S."""

    def Sij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the overlap integral between two STO-nG orbitals."""

        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u in range(n):
            for v in range(m):
                coeff = orb1.cc[u] * orb1.N[u] * orb2.cc[v] * orb2.N[v]
                E_AB = self._E_AB(orb1.gtos[u], orb2.gtos[v])
                prefactor = (math.pi / (orb1.gtos[u].alpha + orb2.gtos[v].alpha)) ** (3 / 2)
                gto1 = orb1.gtos[u]
                gto2 = orb2.gtos[v]
                matrix[u, v] = coeff *E_AB * prefactor * self.integrate_si(gto1, gto2, 'x') * self.integrate_si(gto1, gto2, 'y') * self.integrate_si(gto1, gto2, 'z')

                # coeff = orb1.cc[u] * orb1.N[u] * orb2.cc[v] * orb2.N[v]
                # l1 = orb1.l
    
                # l2 = orb2.l
                
                # if l1 == 0 and l2 == 0:
                #     matrix[u, v] = coeff * self.interate_ss(orb1.gtos[u], orb2.gtos[v])
                # elif l1 == 0 and l2 == 1:
                #     matrix[u, v] = coeff * self.integrate_sp(orb1.gtos[u], orb2.gtos[v])
                # elif l1 == 1 and l2 == 0:
                #     matrix[u, v] = coeff * self.integrate_sp(orb2.gtos[v], orb1.gtos[u])
                # elif l1 == 1 and l2 == 1:
                #     matrix[u, v] = coeff * self.integrate_pp(orb1.gtos[u], orb2.gtos[v])
        return np.sum(matrix)
    
    def _E_AB(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the exponential prefactor E_AB for the overlap integral."""
        prefactor = - gto1.alpha * gto2.alpha / (gto1.alpha + gto2.alpha)
        r_squared = self._r2(gto1, gto2)
        return math.exp(prefactor * r_squared)
    
    def _r2(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the squared distance between the centers of two GTO primitives."""
        dist_x = gto1.coords[0] - gto2.coords[0]
        dist_y = gto1.coords[1] - gto2.coords[1]
        dist_z = gto1.coords[2] - gto2.coords[2]
        return dist_x ** 2 + dist_y ** 2 + dist_z ** 2
    
    def _Pi(self, gto1: STOGPrimitive, gto2: STOGPrimitive, component: Literal['x', 'y', 'z'] = 'x') -> float:
        """Calculate the weighted center and combined exponent for the integral in a given component."""
        alpha_sum = gto1.alpha + gto2.alpha
        if component == 'x':
            return (gto1.alpha * gto1.coords[0] + gto2.alpha * gto2.coords[0]) / alpha_sum
        elif component == 'y':
            return (gto1.alpha * gto1.coords[1] + gto2.alpha * gto2.coords[1]) / alpha_sum
        else:  # component == 'z'
            return (gto1.alpha * gto1.coords[2] + gto2.alpha * gto2.coords[2]) / alpha_sum
    

    def integrate_si(self, gto1: STOGPrimitive, gto2: STOGPrimitive, component: Literal['x', 'y', 'z'] = 'x') -> float:
        """https://content.wolfram.com/sites/19/2012/02/Ho.pdf"""

        A = getattr(gto1, f"coords")[['x', 'y', 'z'].index(component)]
        B = getattr(gto2, f"coords")[['x', 'y', 'z'].index(component)]
        P = self._Pi(gto1, gto2, component)
        
        if getattr(gto1, f"n{component}") == 0 and getattr(gto2, f"n{component}") == 0:
            return 1
        if getattr(gto1, f"n{component}") == 1 and getattr(gto2, f"n{component}") == 0:
            return P - A
        if getattr(gto1, f"n{component}") == 0 and getattr(gto2, f"n{component}") == 1:
            return P - B
        if getattr(gto1, f"n{component}") == 1 and getattr(gto2, f"n{component}") == 1:
            return -(A-P)**2 + (1/(2*(gto1.alpha + gto2.alpha))) + (A-B)*(A-P)


    def interate_ss(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the overlap integral between two s-type GTO primitives."""
        prefactor = math.sqrt(math.pi / (gto1.alpha + gto2.alpha))
        a = gto1.alpha + gto2.alpha
        
        bx = 2*gto1.alpha*gto1.coords[0] + 2*gto2.alpha*gto2.coords[0]
        by = 2*gto1.alpha*gto1.coords[1] + 2*gto2.alpha*gto2.coords[1]
        bz = 2*gto1.alpha*gto1.coords[2] + 2*gto2.alpha*gto2.coords[2]

        cx = gto1.alpha * gto1.coords[0]**2 + gto2.alpha * gto2.coords[0]**2
        cy = gto1.alpha * gto1.coords[1]**2 + gto2.alpha * gto2.coords[1]**2
        cz = gto1.alpha * gto1.coords[2]**2 + gto2.alpha * gto2.coords[2]**2

        resx = math.exp((bx**2 / (4*a)) - cx)
        resy = math.exp((by**2 / (4*a)) - cy)
        resz = math.exp((bz**2 / (4*a)) - cz)
        return prefactor**3 * resx * resy * resz
    
    def integrate_sp(self, gto_s: STOGPrimitive, gto_p: STOGPrimitive) -> float:
        """Calculate the overlap integral between an s-type and a p-type GTO primitive.
        
        Formula: (alpha_s * (X_p - X_s)) / (alpha_s + alpha_p) * integral_ss
        """
        if gto_p.nx == 1:
            center_s = gto_s.coords[0]
            center_p = gto_p.coords[0]
        elif gto_p.ny == 1:
            center_s = gto_s.coords[1]
            center_p = gto_p.coords[1]
        else:
            center_s = gto_s.coords[2]
            center_p = gto_p.coords[2]

        prefactor = (gto_s.alpha * (center_p - center_s)) / (gto_s.alpha + gto_p.alpha)

        return -prefactor * self.interate_ss(gto_s, gto_p)
    
    def integrate_pp(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the overlap integral between two p-type GTO primitives.
        
        For same orbital type (e.g., px-px):
            [1/(2(alpha1+alpha2)) - alpha1*alpha2*(Xi-Xj)^2/(alpha1+alpha2)^2] * integral_ss
        
        For orthogonal orbitals (e.g., px-py):
            [alpha1*(Xi-Xj)/(alpha1+alpha2)] * [alpha1*(Yi-Yj)/(alpha1+alpha2)] * integral_ss
        """
        alpha_sum = gto1.alpha + gto2.alpha
        ss_integral = self.interate_ss(gto1, gto2)
        
        # Check if same orbital type (e.g., px-px, py-py, pz-pz)
        if (gto1.nx == gto2.nx == 1) or (gto1.ny == gto2.ny == 1) or (gto1.nz == gto2.nz == 1):
            # Get the relevant coordinate difference
            if gto1.nx == 1:
                diff = gto1.coords[0] - gto2.coords[0]
            elif gto1.ny == 1:
                diff = gto1.coords[1] - gto2.coords[1]
            else:  # nz == 1
                diff = gto1.coords[2] - gto2.coords[2]
            
            prefactor = (1.0 / (2.0 * alpha_sum)) - (gto1.alpha * gto2.alpha * diff**2) / (alpha_sum**2)
            return prefactor * ss_integral
        
        # Orthogonal orbitals (e.g., px-py, px-pz, py-pz)
        else:
            # Find which coordinates are active
            factor = 1.0
            if gto1.nx == 1:
                diff_x = gto1.coords[0] - gto2.coords[0]
                factor *= (gto1.alpha * diff_x) / alpha_sum
            if gto1.ny == 1:
                diff_y = gto1.coords[1] - gto2.coords[1]
                factor *= (gto1.alpha * diff_y) / alpha_sum
            if gto1.nz == 1:
                diff_z = gto1.coords[2] - gto2.coords[2]
                factor *= (gto1.alpha * diff_z) / alpha_sum
                
            if gto2.nx == 1 and gto1.nx == 0:
                diff_x = gto1.coords[0] - gto2.coords[0]
                factor *= (gto2.alpha * diff_x) / alpha_sum
            if gto2.ny == 1 and gto1.ny == 0:
                diff_y = gto1.coords[1] - gto2.coords[1]
                factor *= (gto2.alpha * diff_y) / alpha_sum
            if gto2.nz == 1 and gto1.nz == 0:
                diff_z = gto1.coords[2] - gto2.coords[2]
                factor *= (gto2.alpha * diff_z) / alpha_sum
            
            return factor * ss_integral


class TijIntegrator:
    """Compute Tij for the kinetic energy matrix T."""

    def __init__(self):
        self.sij = SijIntegrator()

    def Tij(self, orb1: STOGOrbital, orb2: STOGOrbital) -> float:
        """Calculate the kinetic energy integral between two STO-nG orbitals."""

        n = len(orb1.cc)
        m = len(orb2.cc)
        matrix = np.zeros((n, m))
        for u in range(n):
            for v in range(m):
                coeff = orb1.cc[u] * orb1.N[u] * orb2.cc[v] * orb2.N[v]
                l1 = orb1.gtos[u].l
                l2 = orb2.gtos[v].l
                
                if l1 == 0 and l2 == 0:
                    matrix[u, v] = coeff * self.integrate_ss(orb1.gtos[u], orb2.gtos[v])
                elif l1 == 0 and l2 == 1:
                    matrix[u, v] = coeff * self.integrate_sp(orb1.gtos[u], orb2.gtos[v])
                elif l1 == 1 and l2 == 0:
                    matrix[u, v] = coeff * self.integrate_ps(orb1.gtos[u], orb2.gtos[v])
                elif l1 == 1 and l2 == 1:
                    matrix[u, v] = coeff * self.integrate_pp(orb1.gtos[u], orb2.gtos[v])
        return np.sum(matrix)     
    
    def integrate_ss(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the kinetic energy integral between two s-type GTO primitives.
        https://booksite.elsevier.com/9780444594365/downloads/16755_10036.pdf
        """
        prefactor = (gto1.alpha * gto2.alpha) / (gto1.alpha + gto2.alpha)
        R2 = (gto1.coords[0] - gto2.coords[0])**2 + (gto1.coords[1] - gto2.coords[1])**2 + (gto1.coords[2] - gto2.coords[2])**2
        return prefactor * (3 - 2 * prefactor * R2) * self.sij.interate_ss(gto1, gto2)
    
    def integrate_sp(self, gto_s: STOGPrimitive, gto_p: STOGPrimitive) -> float:
        """Calculate the kinetic energy integral between an s-type and a p-type GTO primitive.
        
        Formula: T_sp = alpha_p * S_sp - 2*alpha_p^2 * S_sp(n_p + 2)
        """
        alpha_p = gto_p.alpha
        
        # Base s-p overlap
        s_sp = self.sij.integrate_sp(gto_s, gto_p)
        
        # Create modified p orbital with incremented quantum number
        if gto_p.nx == 1:
            gto_p_mod = STOGPrimitive(coords=gto_p.coords, alpha=gto_p.alpha, nx=3, ny=0, nz=0)
        elif gto_p.ny == 1:
            gto_p_mod = STOGPrimitive(coords=gto_p.coords, alpha=gto_p.alpha, nx=0, ny=3, nz=0)
        else:  # nz == 1
            gto_p_mod = STOGPrimitive(coords=gto_p.coords, alpha=gto_p.alpha, nx=0, ny=0, nz=3)
        
        # Overlap with modified quantum number (need to compute s with p^3 type)
        s_sp_mod = self.sij.integrate_sp(gto_s, gto_p_mod)
        
        return alpha_p * s_sp - 2 * alpha_p**2 * s_sp_mod
    
    def integrate_ps(self, gto_p: STOGPrimitive, gto_s: STOGPrimitive) -> float:
        """Calculate the kinetic energy integral between a p-type and an s-type GTO primitive.
        
        T_ps = -0.5 * integral(p_i * d^2/dr^2(s_j))
        This is different from T_sp due to the non-symmetric nature of the kinetic operator.
        """
        alpha_s = gto_s.alpha
        
        # Base p-s overlap (same as s-p)
        s_ps = self.sij.integrate_sp(gto_s, gto_p)
        
        # For p-s, we differentiate the s orbital which gives simpler result
        # T_ps = 3*alpha_s*S_ps
        return 3 * alpha_s * s_ps
    
    def integrate_pp(self, gto1: STOGPrimitive, gto2: STOGPrimitive) -> float:
        """Calculate the kinetic energy integral between two p-type GTO primitives.
        
        Formula: T_pp = alpha_j*S_pp - 2*alpha_j^2*S_pp(n_j+2) for the active dimension
                        plus contributions from the other two dimensions: alpha_j*S_pp
        """
        alpha_j = gto2.alpha
        
        # Base p-p overlap
        s_pp = self.sij.integrate_pp(gto1, gto2)
        
        # Determine which dimension is active for gto2 and create modified orbital
        if gto2.nx == 1:
            gto2_mod = STOGPrimitive(coords=gto2.coords, alpha=gto2.alpha, nx=3, ny=gto2.ny, nz=gto2.nz)
        elif gto2.ny == 1:
            gto2_mod = STOGPrimitive(coords=gto2.coords, alpha=gto2.alpha, nx=gto2.nx, ny=3, nz=gto2.nz)
        else:  # nz == 1
            gto2_mod = STOGPrimitive(coords=gto2.coords, alpha=gto2.alpha, nx=gto2.nx, ny=gto2.ny, nz=3)
        
        # Overlap with modified quantum number
        s_pp_mod = self.sij.integrate_pp(gto1, gto2_mod)
        
        # For p orbitals: one dimension has angular momentum 1, two have 0
        # T = (alpha_j*3 - 2*alpha_j^2*S_modified) for active dimension
        # Simplification: T_pp = 3*alpha_j*S_pp - 2*alpha_j^2*S_pp(+2)
        return 3 * alpha_j * s_pp - 2 * alpha_j**2 * s_pp_mod
    
class VijIntegrator:
    """Compute Vij for the electron-nuclear attraction matrix V."""

    def Vij(self, orb1: STOGOrbital, orb2: STOGOrbital, atom_coords: tuple[float, float, float], atom_charge: float) -> float:
        """Calculate the electron-nuclear attraction integral between two STO-nG orbitals and a nucleus."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the electron-nuclear attraction integral using the appropriate formulas for STO-nG orbitals.
        return 0.0  # Placeholder return value, which should be replaced with the actual calculation of the electron-nuclear attractio
    
    def integrate_ss(self, gto1: STOGPrimitive, gto2: STOGPrimitive, atom_coords: tuple[float, float, float], atom_charge: float) -> float:
        """Calculate the electron-nuclear attraction integral between two s-type GTO primitives and a nucleus."""
        # This is a placeholder implementation and should be replaced with the actual calculation of the electron-nuclear attraction integral for s-type GTO primitives.
        return 0.0  # Placeholder return value, which should be replaced with the actual calculation of the electron-nuclear attraction integral for s-type GTO primitives.