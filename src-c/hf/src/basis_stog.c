#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../include/basis_stog.h"

double compute_N(double alpha, int *n) {
    /*Calculate the normalization constant for each GTO primitive.*/
    double pi = M_PI;
    double prefactor = pow(2 * alpha / pi, 0.75);
    double numerator = pow(8 * alpha, n[0] + n[1] + n[2]) * tgamma(n[0]+1) * tgamma(n[1]+1) * tgamma(n[2]+1);
    double denominator = tgamma(2*n[0]+1) * tgamma(2*n[1]+1) * tgamma(2*n[2]+1);
    double N = prefactor * sqrt(numerator / denominator);
    return N;
};

/*Function to compute the overlap integral matrix: https://content.wolfram.com/sites/19/2012/02/Ho.pdf*/
double compute_sx(double A, double B, double alpha, double beta, int ai, int bi) {
    /*Compute the integral of the product of two GTOs along a single axis (x, y, or z) using recursion.*/

    double P = (alpha * A + beta * B) / (alpha + beta);

    /*Base cases for recursion*/
    if (ai < 0 || bi < 0) {return 0.0;}
    if (ai == 0 && bi == 0) {return 1;}
    if (ai == 0 && bi == 1) {return -(A-P);}

    /*Index recursion: eq. 11*/
    if (ai > 1 && bi == 0) {
        double term_1 = -(A-P) * compute_sx(A, B, alpha, beta, ai-1, 0);
        double term_2 = (ai-1) / (2 * (alpha + beta)) * compute_sx(A, B, alpha, beta, ai-2, 0);
        return term_1 + term_2;
    }

    /*Transfer equation: eq. 12*/
    double term1 = compute_sx(A, B, alpha, beta, ai+1, bi-1);
    double term2 = (A-B)*compute_sx(A, B, alpha, beta, ai, bi-1);
    return term1 + term2;
}

double compute_Sij(STOOrbital orbital1, STOOrbital orbital2) {
    /*Compute the overlap integral between two STO-nG orbitals by multiplying the integrals along*/

    int n = orbital1.n; 
    int m = orbital2.n;
    double Sij = 0.0;
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOPrimitive gto1 = orbital1.primitives[u];
            STOPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.cc * gto1.N * gto2.cc * gto2.N;
            double r2 = pow(gto1.cords[0] - gto2.cords[0], 2) + pow(gto1.cords[1] - gto2.cords[1], 2) + pow(gto1.cords[2] - gto2.cords[2], 2);
            double E_AB = exp(-gto1.alpha * gto2.alpha * r2 / (gto1.alpha + gto2.alpha));
            double prefactor = pow(M_PI / (gto1.alpha + gto2.alpha), 1.5);
            double sx = compute_sx(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx);
            double sy = compute_sx(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny);
            double sz = compute_sx(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz);
            Sij += coeff * prefactor * E_AB * sx * sy * sz;
        }
    }
    return Sij;
}
