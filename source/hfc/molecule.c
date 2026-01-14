#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include "basis_stog.h"
#include <unistd.h> 

typedef struct  {
    char atom[3]; // Element symbol (e.g., "H", "C", "O")
    int Z; // Atomic number
    double coords[3]; // Atomic coordinates (x, y, z)
    STOOrbital *orbitals; // Array of STO-nG orbitals for the atom
    int num_orbitals; // Number of orbitals for the atom
} Atom;

typedef struct {
    Atom *atoms; // Array of atoms in the molecule
    int num_atoms; // Number of atoms in the molecule
} Molecule;


Molecule parse_molecule_from_file(const char *filename, int num_atoms) {
    /*Parse a molecule from an input file, extracting the atomic symbols, coordinates, and STO-nG orbital information to populate the Molecule struct.*/
    FILE *fptr = fopen(filename, "r");
    if (fptr == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
    }

    char symbol[3];
    double x, y, z;
    int Z;
    int basis;
    Molecule molecule;
    molecule.num_atoms = 0;
    Atom atoms[num_atoms];

    while (fscanf(fptr, "%2s %lf %lf %lf %d %d", symbol, &x, &y, &z, &Z, &basis) == 6) {

        /*Parse atomic information*/
        Atom temp = {
            .atom = *symbol,
            .Z = Z,
            .coords = {x, y, z},
            .orbitals = NULL,
            .num_orbitals = basis,
        };

        /*Parse orbital information*/


        atoms[molecule.num_atoms++] = temp;
    }
    molecule.atoms = atoms;
    molecule.num_atoms = num_atoms;
    return molecule;
}

 
int main(int argc, char *argv[]) {

    int opt;
    char *filename = NULL;
    int num_atoms = 0;

    while((opt = getopt(argc, argv, ":f:n:")) != -1) 
    { 
        switch(opt) 
        { 
            case 'f': 
                printf("Parsing filename: %s\n", optarg); 
                filename = optarg;
                break; 
            case 'n': 
                printf("Number of atoms: %s\n", optarg); 
                num_atoms = atoi(optarg);
                break; 
            case ':': 
                printf("option needs a value\n"); 
                break; 
            case '?': 
                printf("unknown option: %c\n", optopt);
                break; 
        } 
    } 

    Molecule molecule = parse_molecule_from_file(filename, num_atoms);
    printf("Molecule has %d atoms.\n", molecule.num_atoms);
    for (int i = 0; i < molecule.num_atoms; i++) {
        printf("Atom %d: %s, Z=%d, coords=(%.2f, %.2f, %.2f)\n", i+1, molecule.atoms[i].atom, molecule.atoms[i].Z, molecule.atoms[i].coords[0], molecule.atoms[i].coords[1], molecule.atoms[i].coords[2]);
    }
    return 0;
}