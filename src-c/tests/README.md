# Hartree-Fock C Implementation - Unit Tests

This directory contains unit tests for the C implementation of Hartree-Fock calculations.

## Test Files

### test_atom.c
Tests the atom parsing and initialization functionality:
- `test_hydrogen_atom()` - Tests H atom creation
- `test_oxygen_atom()` - Tests O atom with multiple orbitals
- `test_carbon_atom()` - Tests C atom parsing

### test_tensor.c
Tests the 4D tensor implementation for storing 2-electron integrals:
- `test_tensor_allocation()` - Validates memory allocation
- `test_tensor_set_get()` - Tests element access
- `test_tensor_initialization()` - Verifies zero initialization
- `test_tensor_large_dimensions()` - Tests realistic molecular orbital counts
- `test_tensor_symmetry_storage()` - Tests independent storage of permutations

### test_basis.c
Tests basis function integral computations:
- `test_overlap_identical_orbitals()` - Tests <φ|φ> = 1
- `test_kinetic_1s_orbital()` - Tests kinetic energy positivity
- `test_nuclear_attraction()` - Tests attractive potential (negative)
- `test_overlap_symmetry()` - Tests S_ij = S_ji
- `test_kinetic_symmetry()` - Tests T_ij = T_ji

### test_integrals.c
Tests molecular integral matrix computations using H2O molecule from test_data/hf:
- `test_h2o_overlap_matrix()` - Compares computed overlap matrix against reference s.dat
- `test_h2o_kinetic_matrix()` - Compares kinetic energy matrix against reference t.dat
- `test_h2o_nuclear_matrix()` - Compares nuclear attraction matrix against reference v.dat
- `test_h2o_core_hamiltonian()` - Tests H = T + V properties
- `test_h2o_2e_integrals()` - Compares 2-electron integrals against reference eri.dat, tests 8-fold symmetry
- `test_h2o_s_inverse_sqrt()` - Tests S^{-1/2} transformation matrix properties and orthogonality

## Building and Running Tests

### Build all tests:
```bash
make
```

### Run all tests:
```bash
make test
```

### Run individual test suites:
```bash
make test-atom
make test-tensor
make test-basis
make test-integrals
```

### Clean build artifacts:
```bash
make clean
```

## Requirements

- GCC compiler with C11 support
- GSL (GNU Scientific Library) version 2.8 or higher
- Make

## Test Output

Each test prints:
- `PASS: test_name` for successful tests
- `FAIL: test_name (details)` for failed tests
- Final summary: `Tests passed: X/Y`

Exit code:
- 0 if all tests pass
- 1 if any test fails

## Notes

- Tests use a tolerance of 0.1 for comparing against reference data (relaxed for V matrix: 0.5)
- Tests validate physical properties (e.g., kinetic energy > 0, attraction < 0)
- Tests check mathematical properties (symmetry, normalization, orthogonality)
- H2O molecule geometry loaded from `../../test_data/hf/geom.dat`
- Reference data compared from `../../test_data/hf/` directory (s.dat, t.dat, v.dat, eri.dat)
- Helper functions load molecule from geom.dat and matrices from reference files
