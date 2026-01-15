#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include <unistd.h> // For getopt
#include "basis_stog.h"
#include "atom.h"

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
        Atom temp = parse_atom(symbol, x, y, z, basis);
        atoms[molecule.num_atoms++] = temp;
    }
    molecule.atoms = atoms;
    molecule.num_atoms = num_atoms;
    return molecule;
}

 
int main(int argc, char *argv[]) {

    int opt;
    char *filename = NULL;
    char *method = NULL;
    int num_atoms = 0;

    while((opt = getopt(argc, argv, ":f:n:m:")) != -1) 
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
            case 'm': 
                printf("Method: %s\n", optarg); 
                method = optarg;
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

    return 0;
}