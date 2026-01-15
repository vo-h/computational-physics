#ifndef MOLECULE_H // Use include guards to prevent multiple inclusions
#define MOLECULE_H
#include "atom.h"

typedef struct {
    Atom *atoms; // Array of atoms in the molecule
    int num_atoms; // Number of atoms in the molecule
} Molecule;

Molecule parse_molecule_from_file(const char *filename, int num_atoms);
double **compute_S(Molecule *molecule, int num_orbitals);
#endif