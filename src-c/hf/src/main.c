
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h> // For getopt
#include "../include/molecule.h"


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
    printf("Parsed molecule with %d atoms.\n", molecule.num_atoms);
    if (method != NULL && strcmp(method, "overlap") == 0) {
        int num_orbitals = 0;
        for (int i = 0; i < molecule.num_atoms; i++) {
            num_orbitals += molecule.atoms[i].num_orbitals;
        }

        double **S = compute_S(&molecule, num_orbitals);
        printf("Computed overlap matrix S of size %d x %d.\n", num_orbitals, num_orbitals);
        for (int i = 0; i < num_orbitals; i++) {
            for (int j = 0; j < num_orbitals; j++) {
                printf("%f\t", S[i][j]);
            }
            printf("\n");
        }
        // Free the overlap matrix
        for (int i = 0; i < num_orbitals; i++) {
            free(S[i]);
        }
        free(S);
    }

    // Clean up molecule memory
    free_molecule(&molecule);

    return 0;
}