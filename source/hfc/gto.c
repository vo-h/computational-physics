#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <math.h>
#include <stdio.h>

double gto(double *COORS, double *coors, int *n, double alpha) {
    /*Evaluate the Gaussian Type Orbital (GTO) privimitve function at a given point (x, y, z) based on the parameters of the GTO.*/

    // Calculate the normalization constant
    double pi = M_PI;
    double prefactor = pow(2 * alpha / pi, 0.75);
    double numerator = pow(8 * alpha, n[0] + n[1] + n[2]) * tgamma(n[0]+1) * tgamma(n[1]+1) * tgamma(n[2]+1);
    double denominator = tgamma(2*n[0]+1) * tgamma(2*n[1]+1) * tgamma(2*n[2]+1);
    double N = prefactor * sqrt(numerator / denominator);
    
    // Calculate the angular part
    double dist_x = coors[0] - COORS[0];
    double dist_y = coors[1] - COORS[1];
    double dist_z = coors[2] - COORS[2];
    double Y = pow(dist_x, n[0]) * pow(dist_y, n[1]) * pow(dist_z, n[2]);
    
    // Calculate the radial part
    double R = exp(-alpha * (pow(dist_x, 2) + pow(dist_y, 2) + pow(dist_z, 2)));
    return N * Y * R;
}

double sto_ng(double *COORS, double *coors, double* cc, double* alpha, int* nx, int* ny, int* nz, int n) {
    /*Evaluate the STO-nG for an atomic orbital (1s, 2s, 2p, etc.), which is a linear combination of n GTOs.*/
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        int n_vals[3] = {nx[i], ny[i], nz[i]};
        result += cc[i] * gto(COORS, coors, n_vals, alpha[i]);
    }
    return result;
}

