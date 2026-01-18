#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "../include/basis_stog.h"
#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/tensor.h"

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


gsl_matrix *compute_1e_integral(Molecule *molecule, int num_orbitals, char *type) {
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
    gsl_matrix *matrix = gsl_matrix_alloc(num_orbitals, num_orbitals);
    if (matrix == NULL) {
        free(orbitals);
        return NULL;
    }

    // Compute the matrix
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i; j < num_orbitals; j++) {
            if (strcmp(type, "overlap") == 0) {
                gsl_matrix_set(matrix, i, j, compute_Sij(orbitals[i], orbitals[j]));
                gsl_matrix_set(matrix, j, i, gsl_matrix_get(matrix, i, j));
            } else if (strcmp(type, "kinetic") == 0) {
                gsl_matrix_set(matrix, i, j, compute_Tij(orbitals[i], orbitals[j]));
                gsl_matrix_set(matrix, j, i, gsl_matrix_get(matrix, i, j));
            } else if (strcmp(type, "nuclear") == 0) {
                gsl_matrix_set(matrix, i, j, 0.0);
                for (int k = 0; k < molecule->num_atoms; k++) {
                    double R[3] = {molecule->atoms[k].coords[0], molecule->atoms[k].coords[1], molecule->atoms[k].coords[2]};
                    double VijR = compute_VijR(orbitals[i], orbitals[j], R);
                    gsl_matrix_set(matrix, i, j, gsl_matrix_get(matrix, i, j) + -molecule->atoms[k].Z * VijR);
                }
                gsl_matrix_set(matrix, j, i, gsl_matrix_get(matrix, i, j));
            }
        }
    }
    free(orbitals);
    return matrix;
}

tensor4d *compute_2e_integral(Molecule *molecule, int num_orbitals) {
    tensor4d *Vijkl = tensor4d_alloc(num_orbitals, num_orbitals, num_orbitals, num_orbitals);
    /*Compute the 2-electron integral tensor for a set of STO-nG orbitals*/

    /* Get list of orbitals - fixed indexing to use running counter */
    STOOrbital *orbitals = (STOOrbital *)malloc(num_orbitals * sizeof(STOOrbital));
    int orbital_idx = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {
        for (int j = 0; j < molecule->atoms[i].num_orbitals; j++) {
            orbitals[orbital_idx++] = molecule->atoms[i].orbitals[j];
        }
    }

    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            for (int k = 0; k < num_orbitals; k++) {
                for (int l = 0; l < num_orbitals; l++) {
                    if ((i*num_orbitals + j) >= (k*num_orbitals + l)) {
                        double val = compute_Vijkl(orbitals[i], orbitals[j], orbitals[k], orbitals[l]);
                        // 8-fold symmetry: (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (kl|ji) = (lk|ij) = (lk|ji)
                        tensor4d_set(Vijkl, i, j, k, l, val);
                        tensor4d_set(Vijkl, j, i, k, l, val);
                        tensor4d_set(Vijkl, i, j, l, k, val);
                        tensor4d_set(Vijkl, j, i, l, k, val);
                        tensor4d_set(Vijkl, k, l, i, j, val);
                        tensor4d_set(Vijkl, k, l, j, i, val);
                        tensor4d_set(Vijkl, l, k, i, j, val);
                        tensor4d_set(Vijkl, l, k, j, i, val);
                    }
                }
            }
        }
    }
    free(orbitals);
    return Vijkl;
}

/* ------ HARTREE-FOCK ROUTINE ------ */

gsl_matrix *compute_S12(Molecule *molecule, int num_orbitals) {
    gsl_matrix *S = compute_1e_integral(molecule, num_orbitals, "overlap");

    // Allocate space for eigenvalues and eigenvectors
    gsl_vector *eval = gsl_vector_alloc(num_orbitals);
    gsl_matrix *evec = gsl_matrix_alloc(num_orbitals, num_orbitals);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(num_orbitals); // Allocate workspace
    gsl_eigen_symmv(S, eval, evec, w); // Compute eigenvalues and eigenvectors
    gsl_eigen_symmv_free(w); // Free workspace
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    
    /* Compute S^-1/2 */
    gsl_matrix *lambda12 = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix_set_zero(lambda12);
    for (size_t i = 0; i < num_orbitals; i++) {
        gsl_matrix_set(lambda12, i, i, 1.0 / sqrt(gsl_vector_get(eval, i)));
    }
    // gsl_matrix *S_inv_sqrt = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, lambda12, 0.0, S); // S^-1/2 = U * D^-1/2 * U^T
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, S, evec, 0.0, lambda12);
    gsl_matrix_free(S);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return lambda12;
}

gsl_matrix *compute_H(Molecule *molecule, int num_orbitals) {
    /*Compute the electronic Hamiltonian*/
    gsl_matrix *nuc = compute_1e_integral(molecule, num_orbitals, "nuclear");
    gsl_matrix *kin = compute_1e_integral(molecule, num_orbitals, "kinetic");
    gsl_matrix *H = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix_add(H, nuc);
    gsl_matrix_add(H, kin);
    gsl_matrix_free(nuc);
    gsl_matrix_free(kin);
    return H;
}

gsl_matrix *compute_F0(Molecule *molecule, int num_orbitals) {
    /*Compute the initial Fock matrix*/
    gsl_matrix *S12 = compute_S12(molecule, num_orbitals);
    gsl_matrix *H = compute_H(molecule, num_orbitals);
    gsl_matrix *temp = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix *F0 = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, S12, H, 0.0, temp); // F0 = S^-1/2 * H * S^-1/2
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, S12, 0.0, F0);
    gsl_matrix_free(S12);
    gsl_matrix_free(H);
    gsl_matrix_free(temp);
    return F0;
}

gsl_matrix *compute_C0(Molecule *molecule, int num_orbitals) {
    /*Compute the initial coefficient matrix C0 by diagonalizing F0*/
    gsl_matrix *F0 = compute_F0(molecule, num_orbitals);
    gsl_vector *eval = gsl_vector_alloc(num_orbitals);
    gsl_matrix *evec = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix *S12 = compute_S12(molecule, num_orbitals);
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(num_orbitals);
    gsl_eigen_symmv(F0, eval, evec, w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);  // Changed to VAL_ASC to match Python's sort by eigenvalue

    gsl_matrix *C0 = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S12, evec, 0.0, C0);
    
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(F0);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(S12);
    return C0; // C0 is the matrix of eigenvectors
}

gsl_matrix *compute_D0(Molecule *molecule, int num_orbitals) {
    /*Compute the initial density matrix D0 from C0*/
    int occ = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {
        occ += molecule->atoms[i].Z;
    }
    occ /= 2; // Assuming closed-shell

    gsl_matrix *C0 = compute_C0(molecule, num_orbitals);
    gsl_matrix_view C0_sub = gsl_matrix_submatrix (C0, 0, 0, num_orbitals, occ);
    gsl_matrix *sub = &C0_sub.matrix;

    gsl_matrix *D0 = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, sub, sub, 0.0, D0); // D0 = C0_occ @ C0_occ.T (no factor of 2, matches Python)
    gsl_matrix_free(C0);
    return D0;
}

double compute_E0(Molecule *molecule, int num_orbitals) {
    /*Compute the initial total energy E0 using D0 and H*/
    gsl_matrix *H = compute_H(molecule, num_orbitals);
    gsl_matrix *D0 = compute_D0(molecule, num_orbitals);
    double E0 = 0.0;
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            E0 += gsl_matrix_get(D0, i, j) * gsl_matrix_get(H, i, j);
        }
    }
    gsl_matrix_free(H);
    gsl_matrix_free(D0);
    return 2*E0;
}

gsl_matrix *compute_F(gsl_matrix *H, gsl_matrix *D, tensor4d *Vee) {
    /*Compute the Fock matrix F using the formula F[μ,ν] = H[μ,ν] + Σ_λσ D[λ,σ] * [2*(μν|λσ) - (μλ|νσ)]*/
    int num_orbitals = H->size1;
    gsl_matrix *F = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix_memcpy(F, H);

    for (int u = 0; u < num_orbitals; u++) {
        for (int v = 0; v < num_orbitals; v++) {
            for (int l = 0; l < num_orbitals; l++) {
                for (int s = 0; s < num_orbitals; s++) {
                    double val = gsl_matrix_get(D, l, s) * (2*tensor4d_get(Vee, u, v, l, s) - tensor4d_get(Vee, u, l, v, s));
                    gsl_matrix_set(F, u, v, gsl_matrix_get(F, u, v) + val);
                }
            }
        }
    }
    return F;
}

double compute_electronic_energy(gsl_matrix *D, gsl_matrix *H, gsl_matrix *F) {
    /*Compute the total energy E using the formula E = 0.5 * Σ_μν D[μ,ν] * (H[μ,ν] + F[μ,ν])*/
    int num_orbitals = H->size1;
    double E = 0.0;
    for (int u = 0; u < num_orbitals; u++) {
        for (int v = 0; v < num_orbitals; v++) {
            double val = gsl_matrix_get(D, u, v) * (gsl_matrix_get(H, u, v) + gsl_matrix_get(F, u, v));
            E += val;
        }
    }
    return E;
}

void execute_closed_shell_hf(Molecule *molecule, double delta, size_t max_iter) {
    /*Main routine to execute closed-shell Hartree-Fock calculations.*/
    int occ = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {
        occ += molecule->atoms[i].Z;
    }
    occ /= 2; // Assuming closed-shell

    int num_orbitals = 0;
    for (int i = 0; i < molecule->num_atoms; i++) {num_orbitals += molecule->atoms[i].num_orbitals;}

    // Compute integrals once (they don't change during SCF)
    gsl_matrix *H = compute_H(molecule, num_orbitals);
    gsl_matrix *S12 = compute_S12(molecule, num_orbitals);
    tensor4d *Vee = compute_2e_integral(molecule, num_orbitals);
    
    // Initial guess
    gsl_matrix *D = compute_D0(molecule, num_orbitals);
    double E = compute_E0(molecule, num_orbitals);
    
    printf("Iter %2d: E(elec) = %18.12f Hartree\n", 0, E);

    for (int iter = 1; iter <= max_iter; iter++) {
        // Build Fock matrix
        gsl_matrix *F = compute_F(H, D, Vee);

        // Transform Fock matrix to orthogonal basis: F' = S^-1/2 * F * S^-1/2
        gsl_matrix *temp = gsl_matrix_alloc(num_orbitals, num_orbitals);
        gsl_matrix *F_prime = gsl_matrix_alloc(num_orbitals, num_orbitals);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, S12, F, 0.0, temp);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, S12, 0.0, F_prime);
        gsl_matrix_free(temp);

        // Diagonalize F'
        gsl_vector *eval = gsl_vector_alloc(num_orbitals);
        gsl_matrix *evec = gsl_matrix_alloc(num_orbitals, num_orbitals);
        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(num_orbitals);
        gsl_eigen_symmv(F_prime, eval, evec, w);
        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
        gsl_eigen_symmv_free(w);
        gsl_matrix_free(F_prime);

        // Transform eigenvectors back to AO basis: C = S^-1/2 * C'
        gsl_matrix *C = gsl_matrix_alloc(num_orbitals, num_orbitals);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S12, evec, 0.0, C);
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        
        // Build new density matrix from occupied orbitals
        gsl_matrix_view C_occ = gsl_matrix_submatrix(C, 0, 0, num_orbitals, occ);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &C_occ.matrix, &C_occ.matrix, 0.0, D);
        gsl_matrix_free(C);

        // Compute new electronic energy
        double E_new = compute_electronic_energy(D, H, F);
        double delta_E = E_new - E;
        
        printf("Iter %2d: E(elec) = %18.12f Hartree, ΔE = %18.12e\n", iter, E_new, delta_E);

        // Check for convergence
        if (fabs(delta_E) < delta) {
            printf("\nSCF converged after %d iterations!\n", iter);
            printf("Final electronic energy: %.12f Hartree\n", E_new);
            gsl_matrix_free(F);
            break;
        }

        E = E_new;
        gsl_matrix_free(F);
    }

    // Clean up
    gsl_matrix_free(H);
    gsl_matrix_free(S12);
    gsl_matrix_free(D);
    tensor4d_free(Vee);
}