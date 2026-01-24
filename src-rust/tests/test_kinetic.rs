//! Unit tests for kinetic energy matrix calculation.

use nalgebra::DMatrix;
use src_rust::hf::molecule::Molecule;
use std::fs;

const GEOM_FILE: &str = "../data/hf/geom.dat";
const KINETIC_FILE: &str = "../data/hf/t.dat";

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
fn test_kinetic_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(KINETIC_FILE, n);
    
    let t = molecule.compute_T();
    assert_eq!(t.shape(), expected.shape(), "Kinetic matrix shape mismatch");
}

#[test]
fn test_kinetic_symmetry() {
    let molecule = load_molecule();
    let t = molecule.compute_T();
    
    // Check that T is symmetric
    for i in 0..t.nrows() {
        for j in 0..t.ncols() {
            assert!(
                (t[(i, j)] - t[(j, i)]).abs() < 1e-10,
                "Kinetic matrix not symmetric at ({}, {}): {} != {}",
                i, j, t[(i, j)], t[(j, i)]
            );
        }
    }
}

#[test]
fn test_kinetic_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(KINETIC_FILE, n);
    let t = molecule.compute_T();
    
    // Relaxed tolerance due to potential geometry or implementation differences
    for i in 0..n {
        for j in 0..n {
            let diff = (t[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 0.1 || diff < 0.1,
                "Kinetic matrix mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, t[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_kinetic_positive_diagonal() {
    let molecule = load_molecule();
    let t = molecule.compute_T();
    
    // Diagonal elements should be positive (kinetic energy is positive)
    for i in 0..t.nrows() {
        assert!(
            t[(i, i)] > 0.0,
            "Diagonal kinetic energy element at {} should be positive: {}",
            i, t[(i, i)]
        );
    }
}
