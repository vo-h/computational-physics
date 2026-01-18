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
    char basis[10];
    Molecule molecule;
    molecule.num_atoms = 0;
    
    // Allocate atoms array on the heap, not the stack
    Atom *atoms = (Atom *)malloc(num_atoms * sizeof(Atom));
    if (atoms == NULL) {
        fprintf(stderr, "Failed to allocate memory for atoms\n");
        fclose(fptr);
        exit(1);
    }

    while (fscanf(fptr, "%2s %lf %lf %lf %d %s", symbol, &x, &y, &z, &Z, basis) == 6) {
        /*Parse atomic information*/
        if (molecule.num_atoms >= num_atoms) {
            fprintf(stderr, "Warning: File contains more atoms than specified (%d)\n", num_atoms);
            break;
        }
        
        /* Extract basis set information (e.g., "STO-3G") and parse the number of GTO primitives (n) from it. */
        int n;
        printf("Parsing atom: %s at (%.2f, %.2f, %.2f) with basis %s\n", symbol, x, y, z, basis);
        if (sscanf(basis, "STO-%dG", &n) != 1) {
            fprintf(stderr, "Warning: Invalid basis set format for atom %s. Expected 'STO-nG'. Skipping.\n", symbol);
            continue;
        }
        Atom temp = parse_atom(symbol, x, y, z, n);
        atoms[molecule.num_atoms++] = temp;
    }
    
    fclose(fptr);
    molecule.atoms = atoms;
    molecule.num_atoms = num_atoms;
    return molecule;
}

double **compute_1e_integral(Molecule *molecule, int num_orbitals, char *type) {
    /*Compute the full overlap integral matrix for a set of STO-nG orbitals by iterating over all pairs of orbitals and calculating their overlap integrals using compute_Sij().*/
    STOOrbital *orbitals = (STOOrbital *)malloc(num_orbitals * sizeof(STOOrbital));
    
    /* Get list of orbitals - fixed indexing to use running counter */
    int orbital_idx = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {
        for (int j = 0; j < molecule->atoms[i].num_orbitals; j++) {
            orbitals[orbital_idx++] = molecule->atoms[i].orbitals[j];
        }
    }

    /* Allocate matrix for 1-electron integrals */
    double **matrix = (double **)malloc(num_orbitals * sizeof(double *));
    if (matrix == NULL) {
        free(orbitals);
        return NULL;
    }

    for (int i = 0; i < num_orbitals; i++) {
        matrix[i] = (double *)malloc(num_orbitals * sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for 1-electron integral matrix row %d\n", i);
            for (int k = 0; k < i; k++) {
                free(matrix[k]);
            }
            free(matrix);
            free(orbitals);
            return NULL;
        }
    }

    // Initialize the matrix (example)
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i; j < num_orbitals; j++) {
            if (strcmp(type, "overlap") == 0) {
                matrix[i][j] = compute_Sij(orbitals[i], orbitals[j]);
                matrix[j][i] = matrix[i][j];
            } else if (strcmp(type, "kinetic") == 0) {
                matrix[i][j] = compute_Tij(orbitals[i], orbitals[j]);
                matrix[j][i] = matrix[i][j];
            } else if (strcmp(type, "nuclear") == 0) {
                matrix[i][j] = 0.0;
                for (int k = 0; k < molecule->num_atoms; k++) {
                    double R[3] = {molecule->atoms[k].coords[0], molecule->atoms[k].coords[1], molecule->atoms[k].coords[2]};
                    double VijR = compute_VijR(orbitals[i], orbitals[j], R);
                    matrix[i][j] += -molecule->atoms[k].Z * VijR;
                }
                matrix[j][i] = matrix[i][j];
            }
        }
    }
    free(orbitals);
    return matrix;
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
