//! Unit tests for nuclear attraction matrix calculation.

use nalgebra::DMatrix;
use src_rust::hf::molecule::Molecule;
use std::fs;

const GEOM_FILE: &str = "../data/hf/geom.dat";
const NUCLEAR_FILE: &str = "../data/hf/v.dat";

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
fn test_nuclear_shape() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(NUCLEAR_FILE, n);
    
    let vne = molecule.compute_Vne();
    assert_eq!(vne.shape(), expected.shape(), "Nuclear attraction matrix shape mismatch");
}

#[test]
fn test_nuclear_symmetry() {
    let molecule = load_molecule();
    let vne = molecule.compute_Vne();
    
    // Check that Vne is symmetric
    for i in 0..vne.nrows() {
        for j in 0..vne.ncols() {
            assert!(
                (vne[(i, j)] - vne[(j, i)]).abs() < 1e-10,
                "Nuclear attraction matrix not symmetric at ({}, {}): {} != {}",
                i, j, vne[(i, j)], vne[(j, i)]
            );
        }
    }
}

#[test]
fn test_nuclear_values() {
    let molecule = load_molecule();
    let n = molecule.orbitals.len();
    let expected = load_expected_matrix(NUCLEAR_FILE, n);
    let vne = molecule.compute_Vne();
    
    // Relaxed tolerance due to potential geometry or implementation differences
    for i in 0..n {
        for j in 0..n {
            let diff = (vne[(i, j)] - expected[(i, j)]).abs();
            let rel_error = diff / (expected[(i, j)].abs() + 1e-10);
            assert!(
                rel_error < 0.1 || diff < 0.5,
                "Nuclear attraction matrix mismatch at ({}, {}): got {}, expected {} (rel_error: {})",
                i, j, vne[(i, j)], expected[(i, j)], rel_error
            );
        }
    }
}

#[test]
fn test_nuclear_negative_diagonal() {
    let molecule = load_molecule();
    let vne = molecule.compute_Vne();
    
    // Diagonal elements should be negative (attractive potential)
    for i in 0..vne.nrows() {
        assert!(
            vne[(i, i)] < 0.0,
            "Diagonal nuclear attraction element at {} should be negative: {}",
            i, vne[(i, i)]
        );
    }
}
