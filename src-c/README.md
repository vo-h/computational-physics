# Setup

```sh
git clone https://github.com/vo-h/computational-physics.git
cd computational-physics/src-c
```

# Build

To build the Hartree-Fock executable:

```sh
make hfc
```

To build and run tests:

```sh
make test
```

# Usage

After building, the `hf.exe` executable is available in `hf/`:

## Compute 1-electron Matrices

Compute overlap, kinetic energy, and nuclear attraction matrices:

```sh
./hf/hf.exe -f ../data/hf/geom.dat -n 3 -m 1e-mat
```

Output:
- Overlap matrix (S)
- Kinetic energy matrix (T)  
- Nuclear attraction matrix (V)

## View Initial Guesses

Compute the S^(-1/2) orthogonalization matrix:

```sh
./hf/hf.exe -f ../data/hf/geom.dat -n 3 -m guesses
```

## Run Full SCF Calculation

Perform a complete closed-shell Hartree-Fock calculation:

```sh
./hf/hf.exe -f ../data/hf/geom.dat -n 3 -m chf
```

## Command Line Options

- `-f <file>`: Input geometry file (required)
- `-n <number>`: Number of atoms in the molecule (required)
- `-m <method>`: Calculation method (required)
  - `1e-mat`: Compute 1-electron integrals
  - `guesses`: Compute initial guess matrices
  - `chf`: Run closed-shell Hartree-Fock SCF

## Input File Format

Geometry files should contain one atom per line with the format:
```
SYMBOL  X  Y  Z  ATOMIC_NUMBER  BASIS_SET
```

Example (H₂O molecule with 3 atoms):
```
O   0.000000000000  -0.143225816552   0.000000000000   8   STO-3G
H   1.638036840407   1.136548822547  -0.000000000000   1   STO-3G
H  -1.638036840407   1.136548822547  -0.000000000000   1   STO-3G
```