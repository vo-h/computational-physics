#ifndef MOLECULE_H // Use include guards to prevent multiple inclusions
#define MOLECULE_H
#include "atom.h"
#include "../include/tensor.h"
#include <gsl/gsl_matrix.h>

typedef struct {
    Atom *atoms; // Array of atoms in the molecule
    int num_atoms; // Number of atoms in the molecule
} Molecule;
void free_molecule(Molecule *molecule);

Molecule parse_molecule_from_file(const char *filename, int num_atoms);
gsl_matrix *compute_1e_integral(Molecule *molecule, int num_orbitals, char *type);
tensor4d *compute_2e_integral(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_S12(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_H(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_F0(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_C0(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_D0(Molecule *molecule, int num_orbitals);
gsl_matrix *compute_F(gsl_matrix *H, gsl_matrix *D, tensor4d *Vee);
double compute_E0(Molecule *molecule, int num_orbitals);
void execute_closed_shell_hf(Molecule *molecule, double delta, size_t max_iter);

#endif