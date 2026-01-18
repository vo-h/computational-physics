#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../hf/include/tensor.h"

#define TOLERANCE 1e-10
#define ASSERT_NEAR(a, b, msg) \
    if (fabs((a) - (b)) > TOLERANCE) { \
        printf("FAIL: %s (expected %.10f, got %.10f)\n", msg, (double)(b), (double)(a)); \
        return 0; \
    }

#define ASSERT_EQ(a, b, msg) \
    if ((a) != (b)) { \
        printf("FAIL: %s (expected %zu, got %zu)\n", msg, (size_t)(b), (size_t)(a)); \
        return 0; \
    }

int test_tensor_allocation() {
    tensor4d *t = tensor4d_alloc(2, 3, 4, 5);
    
    ASSERT_EQ(t->n1, 2, "Dimension n1");
    ASSERT_EQ(t->n2, 3, "Dimension n2");
    ASSERT_EQ(t->n3, 4, "Dimension n3");
    ASSERT_EQ(t->n4, 5, "Dimension n4");
    
    // Block size should be n1 * n2 * n3 * n4
    size_t expected_size = 2 * 3 * 4 * 5;
    ASSERT_EQ(t->block->size, expected_size, "Block size");
    
    tensor4d_free(t);
    
    printf("PASS: test_tensor_allocation\n");
    return 1;
}

int test_tensor_set_get() {
    tensor4d *t = tensor4d_alloc(3, 3, 3, 3);
    
    // Set some values
    tensor4d_set(t, 0, 0, 0, 0, 1.5);
    tensor4d_set(t, 1, 2, 0, 1, 3.14159);
    tensor4d_set(t, 2, 2, 2, 2, -2.71828);
    
    // Get and verify
    double val1 = tensor4d_get(t, 0, 0, 0, 0);
    double val2 = tensor4d_get(t, 1, 2, 0, 1);
    double val3 = tensor4d_get(t, 2, 2, 2, 2);
    
    ASSERT_NEAR(val1, 1.5, "Get value at (0,0,0,0)");
    ASSERT_NEAR(val2, 3.14159, "Get value at (1,2,0,1)");
    ASSERT_NEAR(val3, -2.71828, "Get value at (2,2,2,2)");
    
    tensor4d_free(t);
    
    printf("PASS: test_tensor_set_get\n");
    return 1;
}

int test_tensor_initialization() {
    tensor4d *t = tensor4d_alloc(2, 2, 2, 2);
    
    // All elements should be initialized to 0
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    double val = tensor4d_get(t, i, j, k, l);
                    ASSERT_NEAR(val, 0.0, "Initial value should be zero");
                }
            }
        }
    }
    
    tensor4d_free(t);
    
    printf("PASS: test_tensor_initialization\n");
    return 1;
}

int test_tensor_large_dimensions() {
    // Test with larger dimensions typical for molecular calculations
    tensor4d *t = tensor4d_alloc(7, 7, 7, 7);
    
    ASSERT_EQ(t->n1, 7, "Large dimension n1");
    ASSERT_EQ(t->block->size, 7 * 7 * 7 * 7, "Large block size");
    
    // Set and get corner values
    tensor4d_set(t, 0, 0, 0, 0, 1.0);
    tensor4d_set(t, 6, 6, 6, 6, 2.0);
    tensor4d_set(t, 3, 3, 3, 3, 3.0);
    
    ASSERT_NEAR(tensor4d_get(t, 0, 0, 0, 0), 1.0, "Corner (0,0,0,0)");
    ASSERT_NEAR(tensor4d_get(t, 6, 6, 6, 6), 2.0, "Corner (6,6,6,6)");
    ASSERT_NEAR(tensor4d_get(t, 3, 3, 3, 3), 3.0, "Middle element");
    
    tensor4d_free(t);
    
    printf("PASS: test_tensor_large_dimensions\n");
    return 1;
}

int test_tensor_symmetry_storage() {
    // Test that different permutations can be set independently
    tensor4d *t = tensor4d_alloc(4, 4, 4, 4);
    
    // Set 8-fold symmetric elements
    tensor4d_set(t, 0, 1, 2, 3, 1.23);
    tensor4d_set(t, 1, 0, 2, 3, 2.34);
    tensor4d_set(t, 0, 1, 3, 2, 3.45);
    tensor4d_set(t, 2, 3, 0, 1, 4.56);
    
    // These should be independent (not automatically symmetric)
    ASSERT_NEAR(tensor4d_get(t, 0, 1, 2, 3), 1.23, "Element (0,1,2,3)");
    ASSERT_NEAR(tensor4d_get(t, 1, 0, 2, 3), 2.34, "Element (1,0,2,3)");
    ASSERT_NEAR(tensor4d_get(t, 0, 1, 3, 2), 3.45, "Element (0,1,3,2)");
    ASSERT_NEAR(tensor4d_get(t, 2, 3, 0, 1), 4.56, "Element (2,3,0,1)");
    
    tensor4d_free(t);
    
    printf("PASS: test_tensor_symmetry_storage\n");
    return 1;
}

int main() {
    printf("Running Tensor Tests...\n");
    printf("======================\n");
    
    int passed = 0;
    int total = 5;
    
    passed += test_tensor_allocation();
    passed += test_tensor_set_get();
    passed += test_tensor_initialization();
    passed += test_tensor_large_dimensions();
    passed += test_tensor_symmetry_storage();
    
    printf("======================\n");
    printf("Tests passed: %d/%d\n", passed, total);
    
    return (passed == total) ? 0 : 1;
}
