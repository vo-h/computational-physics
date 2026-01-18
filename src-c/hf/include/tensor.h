#ifndef TENSOR_H // Use include guards to prevent multiple inclusions
#define TENSOR_H
#include <gsl/gsl_block.h>

// 4D tensor: dimensions [n1][n2][n3][n4]
typedef struct {
    gsl_block *block;
    size_t n1, n2, n3, n4;
} tensor4d;

tensor4d* tensor4d_alloc(size_t n1, size_t n2, size_t n3, size_t n4);
double tensor4d_get(const tensor4d *t, size_t i, size_t j, size_t k, size_t l);
void tensor4d_set(tensor4d *t, size_t i, size_t j, size_t k, size_t l, double val);
void tensor4d_free(tensor4d *t);
#endif