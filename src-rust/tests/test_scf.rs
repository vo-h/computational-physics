//! Unit tests for the SCF (Self-Consistent Field) process.

use nalgebra::DMatrix;
use src_rust::hf::molecule::Molecule;
use std::fs;

const GEOM_FILE: &str = "../data/hf/geom.dat";
const S12_FILE: &str = "../data/hf/scf_s12.dat";
const D0_FILE: &str = "../data/hf/scf_d0.dat";
const F0_FILE: &str = "../data/hf/scf_f0.dat";
const F1_FILE: &str = "../data/hf/scf_f1.dat";

fn load_molecule() -> Molecule {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    Molecule::from_string(&geom)
}

fn load_expected_matrix(filename: &str, n: usize) -> DMatrix<f64> {
    let data = fs::read_to_string(filename)
        .expect("Failed to read data file");
    
    let values: Vec<f64> = data
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    
    DMatrix::from_row_slice(n, n, &values)
}

#[test]
fn test_s12_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(S12_FILE, n);
    
    let s12 = molecule.compute_S12();
    assert_eq!(s12.shape(), expected.shape(), "S^-1/2 matrix shape mismatch");
}

#[test]
fn test_s12_symmetry() {
    let molecule = load_molecule();
    let s12 = molecule.compute_S12();
    
    // Check that S^-1/2 is symmetric
    for i in 0..s12.nrows() {
        for j in 0..s12.ncols() {
            assert!(
                (s12[(i, j)] - s12[(j, i)]).abs() < 1e-10,
                "S^-1/2 matrix not symmetric at ({}, {}): {} != {}",
                i, j, s12[(i, j)], s12[(j, i)]
            );
        }
    }
}

#[test]
fn test_s12_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(S12_FILE, n);
    let s12 = molecule.compute_S12();
    
    // S^-1/2 is computed from eigendecomposition
    for i in 0..n {
        for j in 0..n {
            let diff = (s12[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 1e-6 || diff < 1e-7,
                "S^-1/2 mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, s12[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_s12_orthogonalization() {
    let molecule = load_molecule();
    let s12 = molecule.compute_S12();
    let s = molecule.compute_S();
    
    // S^-1/2^T * S * S^-1/2 = I (identity matrix)
    let result = &s12.transpose() * &s * &s12;
    let n = s12.nrows();
    
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (result[(i, j)] - expected).abs() < 1e-10,
                "S^-1/2 orthogonalization failed at ({}, {}): {} != {}",
                i, j, result[(i, j)], expected
            );
        }
    }
}

#[test]
fn test_f0_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(F0_FILE, n);
    
    let f0 = molecule.compute_F0();
    assert_eq!(f0.shape(), expected.shape(), "F0 matrix shape mismatch");
}

#[test]
fn test_f0_symmetry() {
    let molecule = load_molecule();
    let f0 = molecule.compute_F0();
    
    // Check that F0 is symmetric
    for i in 0..f0.nrows() {
        for j in 0..f0.ncols() {
            assert!(
                (f0[(i, j)] - f0[(j, i)]).abs() < 1e-10,
                "F0 matrix not symmetric at ({}, {}): {} != {}",
                i, j, f0[(i, j)], f0[(j, i)]
            );
        }
    }
}

#[test]
fn test_f0_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(F0_FILE, n);
    let f0 = molecule.compute_F0();
    
    // F0 = S^-1/2^T * H * S^-1/2
    for i in 0..n {
        for j in 0..n {
            let diff = (f0[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 1e-6 || diff < 1e-7,
                "F0 mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, f0[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_c0_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    
    let c0 = molecule.compute_C0();
    assert_eq!(c0.nrows(), n, "C0 rows mismatch");
    assert_eq!(c0.ncols(), n, "C0 cols mismatch");
}

#[test]
fn test_c0_orthonormality() {
    let molecule = load_molecule();
    let c0 = molecule.compute_C0();
    let s = molecule.compute_S();
    
    // C0^T * S * C0 should be identity (orthonormal in overlap metric)
    let result = &c0.transpose() * &s * &c0;
    let n = c0.nrows();
    
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (result[(i, j)] - expected).abs() < 1e-9,
                "C0 orthonormality failed at ({}, {}): {} != {}",
                i, j, result[(i, j)], expected
            );
        }
    }
}

#[test]
fn test_d0_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(D0_FILE, n);
    
    let d0 = molecule.compute_D0();
    assert_eq!(d0.shape(), expected.shape(), "D0 matrix shape mismatch");
}

#[test]
fn test_d0_symmetry() {
    let molecule = load_molecule();
    let d0 = molecule.compute_D0();
    
    // Check that D0 is symmetric
    for i in 0..d0.nrows() {
        for j in 0..d0.ncols() {
            assert!(
                (d0[(i, j)] - d0[(j, i)]).abs() < 1e-10,
                "D0 matrix not symmetric at ({}, {}): {} != {}",
                i, j, d0[(i, j)], d0[(j, i)]
            );
        }
    }
}

#[test]
fn test_d0_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(D0_FILE, n);
    let d0 = molecule.compute_D0();
    
    // D0 = C0_occ * C0_occ^T
    for i in 0..n {
        for j in 0..n {
            let diff = (d0[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 1e-6 || diff < 1e-7,
                "D0 mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, d0[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_d0_idempotency() {
    let molecule = load_molecule();
    let d0 = molecule.compute_D0();
    let s = molecule.compute_S();
    
    // D0 * S * D0 * S ≈ D0 * S (idempotency in overlap metric)
    let ds = &d0 * &s;
    let result = &ds * &ds;
    
    for i in 0..d0.nrows() {
        for j in 0..d0.ncols() {
            assert!(
                (result[(i, j)] - ds[(i, j)]).abs() < 1e-5,
                "D0 idempotency failed at ({}, {}): {} != {}",
                i, j, result[(i, j)], ds[(i, j)]
            );
        }
    }
}

#[test]
fn test_d0_trace() {
    let molecule = load_molecule();
    let d0 = molecule.compute_D0();
    let s = molecule.compute_S();
    
    // Tr(D0 * S) = number of occupied orbitals
    let ds = &d0 * &s;
    let trace: f64 = (0..ds.nrows()).map(|i| ds[(i, i)]).sum();
    
    let n_electrons: u8 = molecule.atoms.iter().map(|a| a.Z).sum();
    let n_occupied = ((n_electrons as f64) / 2.0).ceil();
    
    assert!(
        (trace - n_occupied).abs() < 1e-9,
        "D0 trace mismatch: got {}, expected {}",
        trace, n_occupied
    );
}

#[test]
fn test_scf_fock_structure() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected_f1 = load_expected_matrix(F1_FILE, n);
    
    // Compute F1 manually
    let d = molecule.compute_D0();
    let vee = molecule.compute_Vee();
    let s12 = molecule.compute_S12();
    
    // F = H + G, where G comes from density and two-electron integrals
    let f = molecule.compute_F(&d, &vee);
    let f1_ortho = &s12.transpose() * &f * &s12;
    
    // Test shape
    assert_eq!(f1_ortho.shape(), expected_f1.shape(), "F1 shape mismatch");
    
    // Test symmetry
    for i in 0..n {
        for j in 0..n {
            assert!(
                (f1_ortho[(i, j)] - f1_ortho[(j, i)]).abs() < 1e-10,
                "F1 not symmetric at ({}, {})",
                i, j
            );
        }
    }
    
    // Test magnitude is reasonable (within an order of magnitude)
    let max_f1 = f1_ortho.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let max_expected = expected_f1.iter().map(|x| x.abs()).fold(0.0, f64::max);
    
    assert!(
        max_f1 > 0.1 * max_expected,
        "F1 magnitude too small: {} vs expected {}",
        max_f1, max_expected
    );
    assert!(
        max_f1 < 10.0 * max_expected,
        "F1 magnitude too large: {} vs expected {}",
        max_f1, max_expected
    );
}
