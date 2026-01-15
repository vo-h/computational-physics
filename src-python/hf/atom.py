from pydantic import BaseModel, Field, model_validator
from functools import cached_property
from hf.basis_stog import STOGOrbital
from hf import L
import numpy as np
from typing import Literal, Annotated, Self
import re
import basis_set_exchange as bse
from typing import Annotated
from pydantic import validate_call, Field
from mendeleev import element as Element

class Atom(BaseModel):
    """Atom class, which consists of multiple orbitals."""
    atom: str
    coords: tuple[float, float, float] = (0.0, 0.0, 0.0)  # (X, Y, Z)
    coords_units: Literal['Angstrom', 'Bohr'] = 'Angstrom'
    basis: Annotated[str, Field(pattern=r"STO-[2-6]G")] = "STO-3G"

    @cached_property
    def Z(self) -> int:
        """The atomic number of the atom."""
        return Element(self.atom).atomic_number
    
    @cached_property
    def orbitals(self) -> list[STOGOrbital]:
        """List of all atomic orbitals in the molecule."""
        orbitals = []
        occupations = {"S": 1, "P": 3, "D": 5, "F": 7}

        for orbital in self.get_sto_basis(self.atom, self.basis):
            occ = occupations[orbital['type']]
            for ind in range(occ):
                pattern = r"(?<=STO-)(.*?)(?=G)"
                n = int(re.findall(pattern, self.basis)[0])
                orbitals.append(STOGOrbital(
                    type=orbital['type'],
                    atom=self.atom,
                    coords=self.coords,
                    cc=[orbital['gtos'][i]['coeff'] for i in range(n)],
                    alpha=[orbital['gtos'][i]['alpha'] for i in range(n)],
                    nx=[L[orbital['type']][ind][0] for _ in range(n)],
                    ny=[L[orbital['type']][ind][1] for _ in range(n)],
                    nz=[L[orbital['type']][ind][2] for _ in range(n)]
                ))
        return orbitals
    
    @property
    def cc(self) -> np.ndarray:
        """The contracted coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.cc for orb in self.orbitals])

    @property
    def alpha(self) -> np.ndarray:
        """The exponent coefficients of the GTOs in the STO-3G expansion as a numpy array."""
        return np.stack([orb.alpha for orb in self.orbitals])
    
    @model_validator(mode='after')
    def convert_coords(self) -> Self:
        """Convert the coordinates of the atom to Bohr."""
        if self.coords_units == 'Angstrom':
            self.coords = tuple(coord * 1.8897259886 for coord in self.coords)
            self.coords_units = 'Bohr'
        return self
    
    def __call__(self, x: float, y: float, z: float) -> np.ndarray:
        """Evaluate the orbitals of the atom at a given point (x, y, z)."""
        return np.array([orb(x, y, z) for orb in self.orbitals])


    @staticmethod
    @validate_call
    def get_sto_basis(element: str, basis: Annotated[str, Field(pattern=r"STO-[2-6]G")] = "STO-3G") -> list[dict]:
        """Get the basis set for the given elements in the specified format."""
        info = bse.get_basis(basis, elements=[Element(element).atomic_number], fmt='gaussian94', uncontract_spdf=True)
        info = info.split('\n\n')[-1].strip().split('\n')[1:-1]

        # Group the basis set information by orbital type
        orbitals = []
        curr = ''
        for line in info:
            # Case: new orbital type
            if len(curr) == 0 and not line.startswith(' '):
                curr += line 
            # Case: continuation of current orbital type
            elif len(curr) > 0 and line.startswith(' '):
                curr += '\n' + line
            # Case: new orbital type, but we have already collected one
            elif len(curr) > 0 and not line.startswith(' '):
                orbitals.append(curr)
                curr = line
        if len(curr) > 0:
            orbitals.append(curr)


        # Parse the orbital information into a list of dictionaries
        output = []
        for orbital in orbitals:
            curr = {}
            items = orbital.strip().split('\n')
            curr["type"] = items[0].split()[0]
            curr["gtos"] = []
            for item in items[1:]:
                alpha = item.strip().split()[0]
                cc = item.strip().split()[1]
                curr["gtos"].append({
                    "alpha": float(alpha.strip().lower().replace('d', 'e')),
                    "coeff": float(cc.strip().lower().replace('d', 'e'))
                })
            output.append(curr)
        return output