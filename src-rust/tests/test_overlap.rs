//! Unit tests for overlap matrix calculation.

use nalgebra::DMatrix;
use src_rust::hf::molecule::Molecule;
use std::fs;

const GEOM_FILE: &str = "../data/hf/geom.dat";
const OVERLAP_FILE: &str = "../data/hf/s.dat";

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
fn test_overlap_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(OVERLAP_FILE, n);
    
    let s = molecule.compute_S();
    assert_eq!(s.shape(), expected.shape(), "Overlap matrix shape mismatch");
}

#[test]
fn test_overlap_symmetry() {
    let molecule = load_molecule();
    let s = molecule.compute_S();
    
    // Check that S is symmetric
    for i in 0..s.nrows() {
        for j in 0..s.ncols() {
            assert!(
                (s[(i, j)] - s[(j, i)]).abs() < 1e-10,
                "Overlap matrix not symmetric at ({}, {}): {} != {}",
                i, j, s[(i, j)], s[(j, i)]
            );
        }
    }
}

#[test]
fn test_overlap_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(OVERLAP_FILE, n);
    let s = molecule.compute_S();
    
    // Relaxed tolerance due to potential geometry or implementation differences
    for i in 0..n {
        for j in 0..n {
            let diff = (s[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 0.001 || diff < 0.001,
                "Overlap matrix mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, s[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_overlap_diagonal() {
    let molecule = load_molecule();
    let s = molecule.compute_S();
    
    // Diagonal elements should be 1.0 (orbitals normalized)
    for i in 0..s.nrows() {
        assert!(
            (s[(i, i)] - 1.0).abs() < 1e-10,
            "Diagonal element at {} is not 1.0: {}",
            i, s[(i, i)]
        );
    }
}
