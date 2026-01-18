#include <gsl/gsl_block.h>
#include "../include/tensor.h"

// 4D tensor: dimensions [n1][n2][n3][n4]
tensor4d* tensor4d_alloc(size_t n1, size_t n2, size_t n3, size_t n4) {
    tensor4d *t = malloc(sizeof(tensor4d));
    t->n1 = n1;
    t->n2 = n2;
    t->n3 = n3;
    t->n4 = n4;
    t->block = gsl_block_alloc(n1 * n2 * n3 * n4);
    return t;
}

// Access element [i][j][k][l]
double tensor4d_get(const tensor4d *t, size_t i, size_t j, size_t k, size_t l) {
    size_t idx = i * (t->n2 * t->n3 * t->n4) + 
                 j * (t->n3 * t->n4) + 
                 k * t->n4 + 
                 l;
    return t->block->data[idx];
}

void tensor4d_set(tensor4d *t, size_t i, size_t j, size_t k, size_t l, double val) {
    size_t idx = i * (t->n2 * t->n3 * t->n4) + 
                 j * (t->n3 * t->n4) + 
                 k * t->n4 + 
                 l;
    t->block->data[idx] = val;
}

void tensor4d_free(tensor4d *t) {
    gsl_block_free(t->block);
    free(t);
}