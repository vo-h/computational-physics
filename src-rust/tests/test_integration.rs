//! Integration tests and edge cases for the Hartree-Fock implementation.

use src_rust::hf::atom::Atom;
use src_rust::hf::molecule::Molecule;
use std::fs;

const GEOM_FILE: &str = "../data/hf/geom.dat";

#[test]
fn test_molecule_from_string() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    assert_eq!(molecule.atoms.len(), 3, "H2O should have 3 atoms");
    assert!(molecule.orbitals.len() > 0, "Should have orbitals");
}

#[test]
#[ignore] // Requires network access - ignore in CI
fn test_atom_creation() {
    // Test creating a hydrogen atom
    let h_result = Atom::new("H", (0.0, 0.0, 0.0), 1, "3".to_string());
    
    if h_result.is_err() {
        // Network might be unavailable, skip test
        return;
    }
    
    let h = h_result.unwrap();
    assert_eq!(h.symbol, "H");
    assert_eq!(h.Z, 1);
    assert!(h.orbitals.len() > 0, "H atom should have orbitals");
}

#[test]
fn test_nuclear_repulsion() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let vnn = molecule.compute_Vnn();
    
    // Nuclear repulsion should be positive
    assert!(vnn > 0.0, "Nuclear repulsion should be positive: {}", vnn);
    
    // For H2O with reasonable geometry, Vnn should be around 8-10 Hartree
    assert!(vnn > 5.0 && vnn < 15.0, "Nuclear repulsion out of expected range: {}", vnn);
}

#[test]
fn test_core_hamiltonian() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let h = molecule.compute_H();
    let n = molecule.orbitals.len();
    
    // H should be symmetric
    for i in 0..n {
        for j in 0..n {
            assert!(
                (h[(i, j)] - h[(j, i)]).abs() < 1e-10,
                "Core Hamiltonian not symmetric at ({}, {})",
                i, j
            );
        }
    }
    
    // Diagonal elements should be negative (bound states)
    for i in 0..n {
        assert!(
            h[(i, i)] < 0.0,
            "Core Hamiltonian diagonal should be negative at {}: {}",
            i, h[(i, i)]
        );
    }
}

#[test]
fn test_initial_energy() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let e0 = molecule.compute_E0();
    
    // Electronic energy should be negative for bound systems
    assert!(e0 < 0.0, "Initial electronic energy should be negative: {}", e0);
    
    // For H2O with STO-3G, initial energy should be reasonable
    assert!(e0 > -200.0 && e0 < 0.0, "Initial energy out of expected range: {}", e0);
}

#[test]
fn test_matrix_dimensions_consistency() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let n = molecule.orbitals.len();
    
    // All 2D matrices should be n x n
    let s = molecule.compute_S();
    let t = molecule.compute_T();
    let vne = molecule.compute_Vne();
    let h = molecule.compute_H();
    let s12 = molecule.compute_S12();
    let f0 = molecule.compute_F0();
    let c0 = molecule.compute_C0();
    let d0 = molecule.compute_D0();
    
    assert_eq!(s.shape(), (n, n), "S matrix dimension mismatch");
    assert_eq!(t.shape(), (n, n), "T matrix dimension mismatch");
    assert_eq!(vne.shape(), (n, n), "Vne matrix dimension mismatch");
    assert_eq!(h.shape(), (n, n), "H matrix dimension mismatch");
    assert_eq!(s12.shape(), (n, n), "S12 matrix dimension mismatch");
    assert_eq!(f0.shape(), (n, n), "F0 matrix dimension mismatch");
    assert_eq!(c0.shape(), (n, n), "C0 matrix dimension mismatch");
    assert_eq!(d0.shape(), (n, n), "D0 matrix dimension mismatch");
}

#[test]
fn test_electron_count() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let total_electrons: u8 = molecule.atoms.iter().map(|a| a.Z).sum();
    
    // H2O has 10 electrons (8 from O, 1 each from 2 H)
    assert_eq!(total_electrons, 10, "H2O should have 10 electrons");
}

#[test]
fn test_basis_function_count() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    // STO-3G for H2O: O has 5 functions (1s, 2s, 2px, 2py, 2pz)
    //                  H has 1 function (1s)
    // Total: 5 + 1 + 1 = 7 basis functions
    assert_eq!(molecule.orbitals.len(), 7, "H2O with STO-3G should have 7 basis functions");
}

#[test]
fn test_overlap_positive_semidefinite() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let s = molecule.compute_S();
    
    // Eigenvalues of overlap matrix should all be positive
    let eig = s.symmetric_eigen();
    for eigenval in eig.eigenvalues.iter() {
        assert!(
            *eigenval > 1e-10,
            "Overlap matrix should be positive definite, found eigenvalue: {}",
            eigenval
        );
    }
}

#[test]
fn test_fock_matrix_computation() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    let d = molecule.compute_D0();
    let vee = molecule.compute_Vee();
    let f = molecule.compute_F(&d, &vee);
    
    let n = molecule.orbitals.len();
    
    // Fock matrix should be symmetric
    for i in 0..n {
        for j in 0..n {
            assert!(
                (f[(i, j)] - f[(j, i)]).abs() < 1e-10,
                "Fock matrix not symmetric at ({}, {})",
                i, j
            );
        }
    }
}

#[test]
fn test_orbital_primitive_count() {
    let geom = fs::read_to_string(GEOM_FILE)
        .expect("Failed to read geom.dat");
    let molecule = Molecule::from_string(&geom);
    
    // STO-3G means each orbital has 3 primitives
    for orbital in &molecule.orbitals {
        assert_eq!(
            orbital.gtos.len(),
            3,
            "STO-3G orbitals should have 3 primitives, found {}",
            orbital.gtos.len()
        );
    }
}
