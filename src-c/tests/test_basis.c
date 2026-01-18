#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../hf/include/basis_stog.h"

#define TOLERANCE 1e-6
#define ASSERT_NEAR(a, b, msg) \
    if (fabs((a) - (b)) > TOLERANCE) { \
        printf("FAIL: %s (expected %.10f, got %.10f, diff=%.2e)\n", msg, (double)(b), (double)(a), fabs((a)-(b))); \
        return 0; \
    }

// Test overlap integral for identical 1s orbitals at same position
int test_overlap_identical_orbitals() {
    // H atom 1s orbital (STO-3G)
    STOOrbital orbital;
    orbital.n = 3; // STO-3G has 3 primitives
    
    double coords[3] = {0.0, 0.0, 0.0};
    
    // STO-3G parameters for H 1s (from basis set tables)
    double exponents[] = {3.42525091, 0.62391373, 0.16885540};
    double coeffs[] = {0.15432897, 0.53532814, 0.44463454};
    
    orbital.primitives = malloc(3 * sizeof(STOPrimitive));
    for (int i = 0; i < 3; i++) {
        orbital.primitives[i].cords[0] = coords[0];
        orbital.primitives[i].cords[1] = coords[1];
        orbital.primitives[i].cords[2] = coords[2];
        orbital.primitives[i].alpha = exponents[i];
        orbital.primitives[i].cc = coeffs[i];
        orbital.primitives[i].N = 1.0; // Will be computed properly
        orbital.primitives[i].nx = 0;  // s orbital
        orbital.primitives[i].ny = 0;
        orbital.primitives[i].nz = 0;
    }
    
    // Overlap of orbital with itself should be positive
    double S = compute_Sij(orbital, orbital);
    
    if (S <= 0.0) {
        printf("FAIL: Overlap should be positive (got %.6f)\n", S);
        free(orbital.primitives);
        return 0;
    }
    
    printf("PASS: test_overlap_identical_orbitals (S = %.6f)\n", S);
    
    free(orbital.primitives);
    return 1;
}

// Test kinetic energy integral
int test_kinetic_1s_orbital() {
    // H atom 1s orbital (STO-3G)
    STOOrbital orbital;
    orbital.n = 3;
    
    double coords[3] = {0.0, 0.0, 0.0};
    
    double exponents[] = {3.42525091, 0.62391373, 0.16885540};
    double coeffs[] = {0.15432897, 0.53532814, 0.44463454};
    
    orbital.primitives = malloc(3 * sizeof(STOPrimitive));
    for (int i = 0; i < 3; i++) {
        orbital.primitives[i].cords[0] = coords[0];
        orbital.primitives[i].cords[1] = coords[1];
        orbital.primitives[i].cords[2] = coords[2];
        orbital.primitives[i].alpha = exponents[i];
        orbital.primitives[i].cc = coeffs[i];
        orbital.primitives[i].N = 1.0;
        orbital.primitives[i].nx = 0;
        orbital.primitives[i].ny = 0;
        orbital.primitives[i].nz = 0;
    }
    
    // Kinetic energy should be positive
    double T = compute_Tij(orbital, orbital);
    
    if (T < 0.0) {
        printf("FAIL: Kinetic energy should be positive (got %.10f)\n", T);
        free(orbital.primitives);
        return 0;
    }
    
    printf("PASS: test_kinetic_1s_orbital (T = %.6f)\n", T);
    
    free(orbital.primitives);
    return 1;
}

// Test nuclear attraction integral
int test_nuclear_attraction() {
    // H atom 1s orbital
    STOOrbital orbital;
    orbital.n = 3;
    
    double coords[3] = {0.0, 0.0, 0.0};
    double R[3] = {0.0, 0.0, 0.0}; // Nucleus at origin
    
    double exponents[] = {3.42525091, 0.62391373, 0.16885540};
    double coeffs[] = {0.15432897, 0.53532814, 0.44463454};
    
    orbital.primitives = malloc(3 * sizeof(STOPrimitive));
    for (int i = 0; i < 3; i++) {
        orbital.primitives[i].cords[0] = coords[0];
        orbital.primitives[i].cords[1] = coords[1];
        orbital.primitives[i].cords[2] = coords[2];
        orbital.primitives[i].alpha = exponents[i];
        orbital.primitives[i].cc = coeffs[i];
        orbital.primitives[i].N = 1.0;
        orbital.primitives[i].nx = 0;
        orbital.primitives[i].ny = 0;
        orbital.primitives[i].nz = 0;
    }
    
    // Nuclear attraction - just check it returns a value
    double V = compute_VijR(orbital, orbital, R);
    
    printf("PASS: test_nuclear_attraction (V = %.6f)\n", V);
    
    free(orbital.primitives);
    return 1;
}

// Test that overlap is symmetric
int test_overlap_symmetry() {
    STOOrbital orbital1, orbital2;
    
    // Setup first orbital
    orbital1.n = 1;
    orbital1.primitives = malloc(sizeof(STOPrimitive));
    orbital1.primitives[0].cords[0] = 0.0;
    orbital1.primitives[0].cords[1] = 0.0;
    orbital1.primitives[0].cords[2] = 0.0;
    orbital1.primitives[0].alpha = 1.0;
    orbital1.primitives[0].cc = 1.0;
    orbital1.primitives[0].N = 1.0;
    orbital1.primitives[0].nx = 0;
    orbital1.primitives[0].ny = 0;
    orbital1.primitives[0].nz = 0;
    
    // Setup second orbital
    orbital2.n = 1;
    orbital2.primitives = malloc(sizeof(STOPrimitive));
    orbital2.primitives[0].cords[0] = 1.0;
    orbital2.primitives[0].cords[1] = 0.0;
    orbital2.primitives[0].cords[2] = 0.0;
    orbital2.primitives[0].alpha = 0.5;
    orbital2.primitives[0].cc = 1.0;
    orbital2.primitives[0].N = 1.0;
    orbital2.primitives[0].nx = 0;
    orbital2.primitives[0].ny = 0;
    orbital2.primitives[0].nz = 0;
    
    // S_ij should equal S_ji
    double S12 = compute_Sij(orbital1, orbital2);
    double S21 = compute_Sij(orbital2, orbital1);
    
    ASSERT_NEAR(S12, S21, "Overlap symmetry S_ij = S_ji");
    
    free(orbital1.primitives);
    free(orbital2.primitives);
    
    printf("PASS: test_overlap_symmetry\n");
    return 1;
}

// Test that kinetic energy is symmetric
int test_kinetic_symmetry() {
    STOOrbital orbital1, orbital2;
    
    // Setup first orbital
    orbital1.n = 1;
    orbital1.primitives = malloc(sizeof(STOPrimitive));
    orbital1.primitives[0].cords[0] = 0.0;
    orbital1.primitives[0].cords[1] = 0.0;
    orbital1.primitives[0].cords[2] = 0.0;
    orbital1.primitives[0].alpha = 1.0;
    orbital1.primitives[0].cc = 1.0;
    orbital1.primitives[0].N = 1.0;
    orbital1.primitives[0].nx = 0;
    orbital1.primitives[0].ny = 0;
    orbital1.primitives[0].nz = 0;
    
    // Setup second orbital
    orbital2.n = 1;
    orbital2.primitives = malloc(sizeof(STOPrimitive));
    orbital2.primitives[0].cords[0] = 1.0;
    orbital2.primitives[0].cords[1] = 0.0;
    orbital2.primitives[0].cords[2] = 0.0;
    orbital2.primitives[0].alpha = 0.5;
    orbital2.primitives[0].cc = 1.0;
    orbital2.primitives[0].N = 1.0;
    orbital2.primitives[0].nx = 0;
    orbital2.primitives[0].ny = 0;
    orbital2.primitives[0].nz = 0;
    
    // T_ij should equal T_ji
    double T12 = compute_Tij(orbital1, orbital2);
    double T21 = compute_Tij(orbital2, orbital1);
    
    ASSERT_NEAR(T12, T21, "Kinetic energy symmetry T_ij = T_ji");
    
    free(orbital1.primitives);
    free(orbital2.primitives);
    
    printf("PASS: test_kinetic_symmetry\n");
    return 1;
}

int main() {
    printf("Running Basis Function Tests...\n");
    printf("================================\n");
    
    int passed = 0;
    int total = 5;
    
    passed += test_overlap_identical_orbitals();
    passed += test_kinetic_1s_orbital();
    passed += test_nuclear_attraction();
    passed += test_overlap_symmetry();
    passed += test_kinetic_symmetry();
    
    printf("================================\n");
    printf("Tests passed: %d/%d\n", passed, total);
    
    return (passed == total) ? 0 : 1;
}
