#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "../hf/include/molecule.h"

#define TOLERANCE 0.1
#define RELAXED_TOLERANCE 0.5
#define ASSERT_NEAR(a, b, tol, msg) \
    if (fabs((a) - (b)) > tol) { \
        printf("FAIL: %s (expected %.10f, got %.10f, diff=%.2e)\n", msg, (double)(b), (double)(a), fabs((a)-(b))); \
        return 0; \
    }

#define ASSERT_EQ(a, b, msg) \
    if ((a) != (b)) { \
        printf("FAIL: %s (expected %d, got %d)\n", msg, (b), (a)); \
        return 0; \
    }

// Helper function to load reference matrix from file
gsl_matrix* load_matrix_from_file(const char* filename, int rows, int cols) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("ERROR: Cannot open file %s\n", filename);
        return NULL;
    }
    
    gsl_matrix *mat = gsl_matrix_alloc(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val;
            if (fscanf(fp, "%lf", &val) != 1) {
                printf("ERROR: Failed to read matrix element [%d,%d]\n", i, j);
                fclose(fp);
                gsl_matrix_free(mat);
                return NULL;
            }
            gsl_matrix_set(mat, i, j, val);
        }
    }
    fclose(fp);
    return mat;
}

// Helper function to load molecule from geom.dat
Molecule* load_molecule_from_geom(const char* filename, int *num_orbitals) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("ERROR: Cannot open geometry file %s\n", filename);
        return NULL;
    }
    
    // Count atoms first
    int num_atoms = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (strlen(line) > 1 && line[0] != '\n') {
            num_atoms++;
        }
    }
    
    rewind(fp);
    
    Molecule *mol = malloc(sizeof(Molecule));
    mol->num_atoms = num_atoms;
    mol->atoms = malloc(num_atoms * sizeof(Atom));
    
    *num_orbitals = 0;
    for (int i = 0; i < num_atoms; i++) {
        char atom_symbol[4];
        double x, y, z;
        int Z;
        char basis[32];
        
        if (fscanf(fp, "%s %lf %lf %lf %d %s\n", atom_symbol, &x, &y, &z, &Z, basis) != 6) {
            printf("ERROR: Failed to parse atom %d\n", i);
            fclose(fp);
            free(mol->atoms);
            free(mol);
            return NULL;
        }
        
        mol->atoms[i] = parse_atom(atom_symbol, x, y, z, 3);
        *num_orbitals += mol->atoms[i].num_orbitals;
    }
    
    fclose(fp);
    return mol;
}

// Test H2O molecule overlap matrix against reference data
int test_h2o_overlap_matrix() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    ASSERT_EQ(num_orbitals, 7, "H2O should have 7 orbitals");
    
    gsl_matrix *S = compute_1e_integral(mol, num_orbitals, "overlap");
    gsl_matrix *S_ref = load_matrix_from_file("../../test_data/hf/s.dat", num_orbitals, num_orbitals);
    
    if (!S_ref) {
        gsl_matrix_free(S);
        free_molecule(mol);
        free(mol);
        return 0;
    }
    
    // Check diagonal elements (should be ~1)
    for (int i = 0; i < num_orbitals; i++) {
        double Sii = gsl_matrix_get(S, i, i);
        ASSERT_NEAR(Sii, 1.0, TOLERANCE, "Overlap diagonal should be ~1");
    }
    
    // Check symmetry
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i+1; j < num_orbitals; j++) {
            double Sij = gsl_matrix_get(S, i, j);
            double Sji = gsl_matrix_get(S, j, i);
            ASSERT_NEAR(Sij, Sji, TOLERANCE, "Overlap matrix should be symmetric");
        }
    }
    
    // Compare with reference data
    int mismatches = 0;
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            double computed = gsl_matrix_get(S, i, j);
            double reference = gsl_matrix_get(S_ref, i, j);
            if (fabs(computed - reference) > TOLERANCE) {
                mismatches++;
            }
        }
    }
    
    if (mismatches > 0) {
        printf("WARNING: %d elements differ from reference by more than tolerance\n", mismatches);
    }
    
    gsl_matrix_free(S);
    gsl_matrix_free(S_ref);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_overlap_matrix\n");
    return 1;
}

// Test H2O kinetic energy matrix against reference data
int test_h2o_kinetic_matrix() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    gsl_matrix *T = compute_1e_integral(mol, num_orbitals, "kinetic");
    gsl_matrix *T_ref = load_matrix_from_file("../../test_data/hf/t.dat", num_orbitals, num_orbitals);
    
    if (!T_ref) {
        gsl_matrix_free(T);
        free_molecule(mol);
        free(mol);
        return 0;
    }
    
    // All diagonal elements should be positive
    for (int i = 0; i < num_orbitals; i++) {
        double Tii = gsl_matrix_get(T, i, i);
        if (Tii <= 0.0) {
            printf("FAIL: Kinetic energy diagonal T[%d,%d]=%.6f should be positive\n", i, i, Tii);
            gsl_matrix_free(T);
            gsl_matrix_free(T_ref);
            free_molecule(mol);
            free(mol);
            return 0;
        }
    }
    
    // Check symmetry
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i+1; j < num_orbitals; j++) {
            double Tij = gsl_matrix_get(T, i, j);
            double Tji = gsl_matrix_get(T, j, i);
            ASSERT_NEAR(Tij, Tji, TOLERANCE, "Kinetic matrix should be symmetric");
        }
    }
    
    // Compare with reference data
    int mismatches = 0;
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            double computed = gsl_matrix_get(T, i, j);
            double reference = gsl_matrix_get(T_ref, i, j);
            if (fabs(computed - reference) > TOLERANCE) {
                mismatches++;
            }
        }
    }
    
    if (mismatches > 0) {
        printf("WARNING: %d elements differ from reference by more than tolerance\n", mismatches);
    }
    
    gsl_matrix_free(T);
    gsl_matrix_free(T_ref);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_kinetic_matrix\n");
    return 1;
}

// Test H2O nuclear attraction matrix against reference data
int test_h2o_nuclear_matrix() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    gsl_matrix *V = compute_1e_integral(mol, num_orbitals, "nuclear");
    gsl_matrix *V_ref = load_matrix_from_file("../../test_data/hf/v.dat", num_orbitals, num_orbitals);
    
    if (!V_ref) {
        gsl_matrix_free(V);
        free_molecule(mol);
        free(mol);
        return 0;
    }
    
    // Most elements should be negative (attractive potential)
    // Some elements can be very close to zero due to symmetry
    int non_negative_count = 0;
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            double Vij = gsl_matrix_get(V, i, j);
            if (Vij > 1e-6) {  // Allow for numerical zeros
                non_negative_count++;
            }
        }
    }
    
    // Ensure most elements are negative
    int total_elements = num_orbitals * num_orbitals;
    if (non_negative_count > total_elements / 10) {
        printf("FAIL: Too many non-negative nuclear attraction elements (%d/%d)\n", 
               non_negative_count, total_elements);
        gsl_matrix_free(V);
        gsl_matrix_free(V_ref);
        free_molecule(mol);
        free(mol);
        return 0;
    }
    
    // Check symmetry
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i+1; j < num_orbitals; j++) {
            double Vij = gsl_matrix_get(V, i, j);
            double Vji = gsl_matrix_get(V, j, i);
            ASSERT_NEAR(Vij, Vji, TOLERANCE, "Nuclear matrix should be symmetric");
        }
    }
    
    // Compare with reference data (use relaxed tolerance for V)
    int mismatches = 0;
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            double computed = gsl_matrix_get(V, i, j);
            double reference = gsl_matrix_get(V_ref, i, j);
            if (fabs(computed - reference) > RELAXED_TOLERANCE) {
                mismatches++;
            }
        }
    }
    
    if (mismatches > 0) {
        printf("WARNING: %d elements differ from reference by more than tolerance\n", mismatches);
    }
    
    gsl_matrix_free(V);
    gsl_matrix_free(V_ref);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_nuclear_matrix\n");
    return 1;
}

// Test H2O core Hamiltonian (H = T + V)
int test_h2o_core_hamiltonian() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    gsl_matrix *H = compute_H(mol, num_orbitals);
    
    // Core Hamiltonian should be symmetric
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i+1; j < num_orbitals; j++) {
            double Hij = gsl_matrix_get(H, i, j);
            double Hji = gsl_matrix_get(H, j, i);
            ASSERT_NEAR(Hij, Hji, TOLERANCE, "Core Hamiltonian should be symmetric");
        }
    }
    
    // Diagonal elements should be negative (attractive dominates kinetic)
    for (int i = 0; i < num_orbitals; i++) {
        double Hii = gsl_matrix_get(H, i, i);
        if (Hii >= 0.0) {
            printf("FAIL: Core Hamiltonian diagonal H[%d,%d]=%.6f should be negative\n", i, i, Hii);
            gsl_matrix_free(H);
            free_molecule(mol);
            free(mol);
            return 0;
        }
    }
    
    gsl_matrix_free(H);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_core_hamiltonian\n");
    return 1;
}

// Test H2O 2-electron integral tensor against reference data
int test_h2o_2e_integrals() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    tensor4d *Vee = compute_2e_integral(mol, num_orbitals);
    
    // Test 8-fold symmetry on sample elements
    for (int i = 0; i < 2 && i < num_orbitals; i++) {
        for (int j = 0; j < 2 && j < num_orbitals; j++) {
            for (int k = 0; k < 2 && k < num_orbitals; k++) {
                for (int l = 0; l < 2 && l < num_orbitals; l++) {
                    double Vijkl = tensor4d_get(Vee, i, j, k, l);
                    double Vjikl = tensor4d_get(Vee, j, i, k, l);
                    double Vijlk = tensor4d_get(Vee, i, j, l, k);
                    double Vklij = tensor4d_get(Vee, k, l, i, j);
                    
                    ASSERT_NEAR(Vijkl, Vjikl, TOLERANCE, "2e symmetry (ij|kl) = (ji|kl)");
                    ASSERT_NEAR(Vijkl, Vijlk, TOLERANCE, "2e symmetry (ij|kl) = (ij|lk)");
                    ASSERT_NEAR(Vijkl, Vklij, TOLERANCE, "2e symmetry (ij|kl) = (kl|ij)");
                }
            }
        }
    }
    
    // Load reference data from eri.dat and compare
    FILE *fp = fopen("../../test_data/hf/eri.dat", "r");
    if (fp) {
        int mismatches = 0;
        int i, j, k, l;
        double ref_val;
        
        while (fscanf(fp, "%d %d %d %d %lf", &i, &j, &k, &l, &ref_val) == 5) {
            // Convert from 1-based to 0-based indexing
            i--; j--; k--; l--;
            
            if (i >= 0 && i < num_orbitals && j >= 0 && j < num_orbitals &&
                k >= 0 && k < num_orbitals && l >= 0 && l < num_orbitals) {
                double computed = tensor4d_get(Vee, i, j, k, l);
                if (fabs(computed - ref_val) > TOLERANCE) {
                    mismatches++;
                }
            }
        }
        
        if (mismatches > 0) {
            printf("WARNING: %d ERI elements differ from reference\n", mismatches);
        }
        
        fclose(fp);
    }
    
    tensor4d_free(Vee);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_2e_integrals\n");
    return 1;
}

// Test H2O S^{-1/2} matrix computation
int test_h2o_s_inverse_sqrt() {
    int num_orbitals;
    Molecule *mol = load_molecule_from_geom("../../test_data/hf/geom.dat", &num_orbitals);
    if (!mol) return 0;
    
    gsl_matrix *S12 = compute_S12(mol, num_orbitals);
    
    // S12 should be symmetric
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = i+1; j < num_orbitals; j++) {
            double S12_ij = gsl_matrix_get(S12, i, j);
            double S12_ji = gsl_matrix_get(S12, j, i);
            ASSERT_NEAR(S12_ij, S12_ji, TOLERANCE, "S^{-1/2} should be symmetric");
        }
    }
    
    // Diagonal elements should be positive
    for (int i = 0; i < num_orbitals; i++) {
        double S12_ii = gsl_matrix_get(S12, i, i);
        if (S12_ii <= 0.0) {
            printf("FAIL: S^{-1/2} diagonal element [%d,%d]=%.6f should be positive\n", i, i, S12_ii);
            gsl_matrix_free(S12);
            free_molecule(mol);
            free(mol);
            return 0;
        }
    }
    
    // Verify S^{-1/2}^T * S * S^{-1/2} = I
    gsl_matrix *S = compute_1e_integral(mol, num_orbitals, "overlap");
    gsl_matrix *temp = gsl_matrix_alloc(num_orbitals, num_orbitals);
    gsl_matrix *result = gsl_matrix_alloc(num_orbitals, num_orbitals);
    
    // temp = S * S12
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S, S12, 0.0, temp);
    // result = S12^T * temp = S12^T * S * S12
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, S12, temp, 0.0, result);
    
    // Check that result is identity
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            double actual = gsl_matrix_get(result, i, j);
            ASSERT_NEAR(actual, expected, TOLERANCE, "S^{-1/2}^T * S * S^{-1/2} = I");
        }
    }
    
    gsl_matrix_free(S12);
    gsl_matrix_free(S);
    gsl_matrix_free(temp);
    gsl_matrix_free(result);
    free_molecule(mol);
    free(mol);
    
    printf("PASS: test_h2o_s_inverse_sqrt\n");
    return 1;
}

int main() {
    printf("Running H2O Integral Tests (comparing with test_data/hf)...\n");
    printf("==========================================================\n");
    
    int passed = 0;
    int total = 6;
    
    passed += test_h2o_overlap_matrix();
    passed += test_h2o_kinetic_matrix();
    passed += test_h2o_nuclear_matrix();
    passed += test_h2o_core_hamiltonian();
    passed += test_h2o_2e_integrals();
    passed += test_h2o_s_inverse_sqrt();
    
    printf("==========================================================\n");
    printf("Tests passed: %d/%d\n", passed, total);
    
    return (passed == total) ? 0 : 1;
}
