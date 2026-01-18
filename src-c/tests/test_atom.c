#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../hf/include/atom.h"

#define TOLERANCE 1e-6
#define ASSERT_NEAR(a, b, msg) \
    if (fabs((a) - (b)) > TOLERANCE) { \
        printf("FAIL: %s (expected %.10f, got %.10f)\n", msg, (double)(b), (double)(a)); \
        return 0; \
    }

#define ASSERT_EQ(a, b, msg) \
    if ((a) != (b)) { \
        printf("FAIL: %s (expected %d, got %d)\n", msg, (b), (a)); \
        return 0; \
    }

#define ASSERT_STR_EQ(a, b, msg) \
    if (strcmp(a, b) != 0) { \
        printf("FAIL: %s (expected %s, got %s)\n", msg, (b), (a)); \
        return 0; \
    }

int test_hydrogen_atom() {
    Atom h = parse_atom("H", 0.0, 0.0, 0.0, 3);
    
    ASSERT_STR_EQ(h.atom, "H", "Hydrogen symbol");
    ASSERT_EQ(h.Z, 1, "Hydrogen atomic number");
    ASSERT_NEAR(h.coords[0], 0.0, "H x-coordinate");
    ASSERT_NEAR(h.coords[1], 0.0, "H y-coordinate");
    ASSERT_NEAR(h.coords[2], 0.0, "H z-coordinate");
    ASSERT_EQ(h.num_gtos, 3, "Number of GTOs");
    ASSERT_EQ(h.num_orbitals, 1, "Number of orbitals for H");
    
    // Clean up
    for (int i = 0; i < h.num_orbitals; i++) {
        free(h.orbitals[i].primitives);
    }
    free(h.orbitals);
    
    printf("PASS: test_hydrogen_atom\n");
    return 1;
}

int test_oxygen_atom() {
    Atom o = parse_atom("O", 1.0, 2.0, 3.0, 3);
    
    ASSERT_STR_EQ(o.atom, "O", "Oxygen symbol");
    ASSERT_EQ(o.Z, 8, "Oxygen atomic number");
    ASSERT_NEAR(o.coords[0], 1.0, "O x-coordinate");
    ASSERT_NEAR(o.coords[1], 2.0, "O y-coordinate");
    ASSERT_NEAR(o.coords[2], 3.0, "O z-coordinate");
    ASSERT_EQ(o.num_gtos, 3, "Number of GTOs");
    ASSERT_EQ(o.num_orbitals, 5, "Number of orbitals for O (1s, 2s, 2px, 2py, 2pz)");
    
    // Clean up
    for (int i = 0; i < o.num_orbitals; i++) {
        free(o.orbitals[i].primitives);
    }
    free(o.orbitals);
    
    printf("PASS: test_oxygen_atom\n");
    return 1;
}

int test_carbon_atom() {
    Atom c = parse_atom("C", -1.5, 0.5, -0.5, 3);
    
    ASSERT_STR_EQ(c.atom, "C", "Carbon symbol");
    ASSERT_EQ(c.Z, 6, "Carbon atomic number");
    ASSERT_NEAR(c.coords[0], -1.5, "C x-coordinate");
    ASSERT_NEAR(c.coords[1], 0.5, "C y-coordinate");
    ASSERT_NEAR(c.coords[2], -0.5, "C z-coordinate");
    ASSERT_EQ(c.num_orbitals, 5, "Number of orbitals for C");
    
    // Clean up
    for (int i = 0; i < c.num_orbitals; i++) {
        free(c.orbitals[i].primitives);
    }
    free(c.orbitals);
    
    printf("PASS: test_carbon_atom\n");
    return 1;
}

int main() {
    printf("Running Atom Tests...\n");
    printf("====================\n");
    
    int passed = 0;
    int total = 3;
    
    passed += test_hydrogen_atom();
    passed += test_oxygen_atom();
    passed += test_carbon_atom();
    
    printf("====================\n");
    printf("Tests passed: %d/%d\n", passed, total);
    
    return (passed == total) ? 0 : 1;
}
