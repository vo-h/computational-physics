
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h> // For getopt
#include <gsl/gsl_matrix.h>
#include "../include/molecule.h"
#include "../include/tensor.h"

void print_2D_tensor(gsl_matrix *matrix, int num_orbitals) {
    /*Helper function to print a matrix in a readable format.*/
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            printf("%10.6f\t", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
    printf("\n");
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
    printf("Parsed molecule with %d atoms.\n", molecule.num_atoms);
    if (method != NULL && (strcmp(method, "1e-mat") == 0)) {
        
        /* Compute the number of orbitals in the molecule */
        int num_orbitals = 0;
        for (int i = 0; i < molecule.num_atoms; i++) {num_orbitals += molecule.atoms[i].num_orbitals;}

        char *e1_methods[] = {"overlap", "kinetic", "nuclear"};
        for (int i = 0; i < 3; i++) {
            gsl_matrix *matrix = compute_1e_integral(&molecule, num_orbitals, e1_methods[i]);
            if (matrix == NULL) {
                fprintf(stderr, "Error: Failed to compute %s matrix\n",  e1_methods[i]);
                free_molecule(&molecule);
                return 1;
            }
            printf("Computed %s matrix of size %d x %d.\n", e1_methods[i], num_orbitals, num_orbitals);
            print_2D_tensor(matrix, num_orbitals);
            gsl_matrix_free(matrix); // Free the matrix
        }
    } else if (method != NULL && strcmp(method, "guesses") == 0) {
        /* Compute the number of orbitals in the molecule */
        int num_orbitals = 0;
        for (int i = 0; i < molecule.num_atoms; i++) {num_orbitals += molecule.atoms[i].num_orbitals;}

        gsl_matrix *S12 = compute_S12(&molecule, num_orbitals);
        if (S12 == NULL) {
            fprintf(stderr, "Error: Failed to compute S^-1/2 matrix\n");
            free_molecule(&molecule);
            return 1;
        }
        printf("Computed S^-1/2 matrix of size %d x %d.\n", num_orbitals, num_orbitals);
        print_2D_tensor(S12, num_orbitals);
        gsl_matrix_free(S12); // Free the tensor
    } else if (method != NULL && strcmp(method, "chf") == 0) {
        execute_closed_shell_hf(&molecule, 1e-6, 100);
        free_molecule(&molecule);
        return 1;
    }
    
    else{
        fprintf(stderr, "Error: Invalid method specified. Use '1e-mat' or 'guesses'.\n");
        free_molecule(&molecule);
        return 1;
    }

    free_molecule(&molecule); // Clean up molecule memory
    return 0;
}
