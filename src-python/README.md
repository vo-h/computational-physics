# Setup

```sh
git clone https://github.com/vo-h/computational-physics.git
cd computational-physics/src-python
pip install -e .
```

# Usage: Hartree-Fock Calculations

```python
from hf.molecule import Molecule
import numpy as np

# Load molecule from geometry file
molecule = Molecule.from_file("../data/hf/geom.dat")

# Compute 1-electron integrals
S = molecule.S  # Overlap matrix
T = molecule.T  # Kinetic energy matrix
V = molecule.Vne  # Nuclear attraction matrix

print(f"Overlap matrix:\n{S}\n")
print(f"Kinetic energy matrix:\n{T}\n")
print(f"Nuclear attraction matrix:\n{V}\n")

# Perform SCF calculation
molecule.CHF()
print(f"SCF converged. Final energy: {molecule.E:.6f} Hartree")
```

## Command Line Interface

After installation, the `hf` command is available:

```sh
# Compute 1-electron matrices (overlap, kinetic, nuclear attraction)
hf --input-file ../data/hf/geom.dat --method 1e-mat

# View initial guesses (Fock matrix, coefficients, density matrix, energy)
hf --input-file ../data/hf/geom.dat --method guesses

# Run full SCF calculation (default)
hf --input-file ../data/hf/geom.dat --method scf
# or simply:
hf --input-file ../data/hf/geom.dat
```

### Input File Format

Geometry files should contain one atom per line with the format:
```
SYMBOL  X  Y  Z  ATOMIC_NUMBER  BASIS_SET
```

Example (H₂O molecule):
```
O   0.000000000000  -0.143225816552   0.000000000000   8   STO-3G
H   1.638036840407   1.136548822547  -0.000000000000   1   STO-3G
H  -1.638036840407   1.136548822547  -0.000000000000   1   STO-3G
```