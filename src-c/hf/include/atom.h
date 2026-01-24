#ifndef ATOM_H // Use include guards to prevent multiple inclusions
#define ATOM_H
#include "basis_stog.h"

typedef struct  {
    char atom[3]; // Element symbol (e.g., "H", "C", "O")
    int Z; // Atomic number
    double coords[3]; // Atomic coordinates (x, y, z)
    STOGOrbital *orbitals; // Array of STO-nG orbitals for the atom
    int num_orbitals; // Number of orbitals for the atom
    int num_gtos; // Number of GTO primitives per orbital
} Atom;

Atom parse_atom(char *symbol, double x, double y, double z, int num_gtos);
#endif
