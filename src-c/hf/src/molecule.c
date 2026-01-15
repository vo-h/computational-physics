#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../include/basis_stog.h"
#include "../include/atom.h"
#include "../include/molecule.h"

Molecule parse_molecule_from_file(const char *filename, int num_atoms) {
    /*Parse a molecule from an input file, extracting the atomic symbols, coordinates, and STO-nG orbital information to populate the Molecule struct.*/
    FILE *fptr = fopen(filename, "r");
    if (fptr == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(1);
    }

    char symbol[3];
    double x, y, z;
    int Z;
    int basis;
    Molecule molecule;
    molecule.num_atoms = 0;
    
    // Allocate atoms array on the heap, not the stack
    Atom *atoms = (Atom *)malloc(num_atoms * sizeof(Atom));
    if (atoms == NULL) {
        fprintf(stderr, "Failed to allocate memory for atoms\n");
        fclose(fptr);
        exit(1);
    }

    while (fscanf(fptr, "%2s %lf %lf %lf %d %d", symbol, &x, &y, &z, &Z, &basis) == 6) {
        /*Parse atomic information*/
        if (molecule.num_atoms >= num_atoms) {
            fprintf(stderr, "Warning: File contains more atoms than specified (%d)\n", num_atoms);
            break;
        }
        Atom temp = parse_atom(symbol, x, y, z, basis);
        atoms[molecule.num_atoms++] = temp;
    }
    
    fclose(fptr);
    molecule.atoms = atoms;
    molecule.num_atoms = num_atoms;
    return molecule;
}

double **compute_S(Molecule *molecule, int num_orbitals) {
    /*Compute the full overlap integral matrix for a set of STO-nG orbitals by iterating over all pairs of orbitals and calculating their overlap integrals using compute_Sij().*/
    STOOrbital *orbitals = (STOOrbital *)malloc(num_orbitals * sizeof(STOOrbital));
    
    /* Get list of orbitals - fixed indexing to use running counter */
    int orbital_idx = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {
        for (int j = 0; j < molecule->atoms[i].num_orbitals; j++) {
            orbitals[orbital_idx++] = molecule->atoms[i].orbitals[j];
        }
    }

    /* Allocate matrix for overlap integrals */
    double **S = (double **)malloc(num_orbitals * sizeof(double *));
    if (S == NULL) {
        free(orbitals);
        return NULL;
    }

    for (int i = 0; i < num_orbitals; i++) {
        S[i] = (double *)malloc(num_orbitals * sizeof(double));
        if (S[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for overlap matrix row %d\n", i);
            for (int k = 0; k < i; k++) {
                free(S[k]);
            }
            free(S);
            free(orbitals);
            return NULL;
        }
    }

    // Initialize the matrix (example)
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            S[i][j] = compute_Sij(orbitals[i], orbitals[j]);
        }
    }
    free(orbitals);
    return S;
}

void free_molecule(Molecule *molecule) {
    /*Free all memory allocated for the molecule*/
    if (molecule->atoms != NULL) {
        for (int i = 0; i < molecule->num_atoms; i++) {
            if (molecule->atoms[i].orbitals != NULL) {
                for (int j = 0; j < molecule->atoms[i].num_orbitals; j++) {
                    free(molecule->atoms[i].orbitals[j].primitives);
                }
                free(molecule->atoms[i].orbitals);
            }
        }
        free(molecule->atoms);
        molecule->atoms = NULL;
    }
}


