#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_sf_hyperg.h>
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

double compute_E(double A, double B, double alpha, double beta, int ai, int bi, int t) {
    /*Compute the nuclear attraction integral between two GTOs along a single axis using recursion.
      Arguments:
          A, B: origin of Gaussian 'a' and 'b'
          alpha, beta: exponent of Gaussian 'a' and 'b'
          ai, bi: orbital angular momentum number on Gaussian 'a' and 'b'
          t: number nodes in Hermite (depends on type of integral, e.g. always zero for overlap integrals)
    */
    double p = alpha + beta;
    double q = alpha * beta / p;
    double Qx = A - B;

    if (t < 0 || t > ai+bi) {return 0.0;}
    if (ai == 0 && bi == 0 && t == 0) {return exp(-q * Qx * Qx);}
    if (bi == 0) {
        double term1 = (1/(2*p)) * compute_E(A, B, alpha, beta, ai-1, bi, t-1);
        double term2 = (q * Qx / alpha) * compute_E(A, B, alpha, beta, ai-1, bi, t);
        double term3 = (t+1) * compute_E(A, B, alpha, beta, ai-1, bi, t+1);
        return term1 - term2 + term3;
    }

    double term1 = (1/(2*p)) * compute_E(A, B, alpha, beta, ai, bi-1, t-1);
    double term2 = (q * Qx / beta) * compute_E(A, B, alpha, beta, ai, bi-1, t);
    double term3 = (t+1) * compute_E(A, B, alpha, beta, ai, bi-1, t+1);
    return term1 + term2 + term3;
}

double compute_boys(int n, double T) {
    /*Compute the Boys function using the confluent hypergeometric function 1F1, which is used in the evaluation of nuclear attraction integrals.*/
    return gsl_sf_hyperg_1F1(n+0.5, n+1.5, -T) / (2.0*n + 1.0);
}

double compute_R(int t, int u, int v, int n, double p, double P[3], double C[3]) {
    /*Returns the Coulomb auxiliary Hermite integrals
      Arguments:
          t,u,v: order of Coulomb Hermite derivative in x,y,z
          n: order of Boys function
          P: Gaussian composite center P
          C: Nuclear center C
    */
    double RPC = sqrt(pow(P[0]-C[0], 2) + pow(P[1]-C[1], 2) + pow(P[2]-C[2], 2));
    double val = 0.0;
    
    if (t < 0 || u < 0 || v < 0) {return 0.0;}
    if (t == 0 && u == 0 && v == 0) {
        val += pow(-2*p, n) * compute_boys(n, p * RPC * RPC);
    } else if (t == 0 && u == 0) {
        val += (v-1) * compute_R(t, u, v-2, n+1, p, P, C);
        val += (P[2]-C[2]) * compute_R(t, u, v-1, n+1, p, P, C);
    } else if (t == 0) {
        val += (u-1) * compute_R(t, u-2, v, n+1, p, P, C);
        val += (P[1]-C[1]) * compute_R(t, u-1, v, n+1, p, P, C);
    } else {
        val += (t-1) * compute_R(t-2, u, v, n+1, p, P, C);
        val += (P[0]-C[0]) * compute_R(t-1, u, v, n+1, p, P, C);
    }
    return val;
}

double compute_D(double A, double B, double alpha, double beta, int ai, int bi, int t) {
    /*Returns the Dawson function D_n(T) using recursion.
      Arguments:
          A, B: origin of Gaussian 'a' and 'b'
          alpha, beta: exponent of Gaussian 'a' and 'b'
          ai, bi: orbital angular momentum number on Gaussian 'a' and 'b'
    */
    double p = alpha + beta;
    if (t == 0) {return compute_E(A, B, alpha, beta, ai, bi, 0) * pow(M_PI/p, 0.5);}
    double term1 = bi * compute_D(A, B, alpha, beta, ai, bi-1, t-1);
    double term2 = 2 * beta * compute_D(A, B, alpha, beta, ai, bi+1, t-1);
    return term1 - term2;
}

double compute_Sij(STOGOrbital orbital1, STOGOrbital orbital2) {
    /*Compute the overlap integral between two STO-nG orbitals by multiplying the integrals along*/

    int n = orbital1.n; 
    int m = orbital2.n;
    double Sij = 0.0;
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOGPrimitive gto1 = orbital1.primitives[u];
            STOGPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.cc * gto1.N * gto2.cc * gto2.N;
            double prefactor = pow(M_PI / (gto1.alpha + gto2.alpha), 1.5);
            double s_x = compute_E(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, 0); 
            double s_y = compute_E(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, 0);
            double s_z = compute_E(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, 0);
            Sij += coeff * prefactor * s_x * s_y * s_z;
        }
    }
    return Sij;
}

double compute_Tij(STOGOrbital orbital1, STOGOrbital orbital2) {
    /*Compute the kinetic energy integral between two STO-nG orbitals by multiplying the integrals along each axis and summing contributions from all GTO primitives.*/

    int n = orbital1.n; 
    int m = orbital2.n;
    double Tij = 0.0;
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOGPrimitive gto1 = orbital1.primitives[u];
            STOGPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.cc * gto1.N * gto2.cc * gto2.N;
            double term1 = compute_D(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, 2);
            double term2 = compute_D(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, 2);
            double term3 = compute_D(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, 2);
            double term4 = compute_D(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, 0);
            double term5 = compute_D(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, 0);
            double term6 = compute_D(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, 0);
            Tij += -0.5 * coeff * (term1*term5*term6 + term4*term2*term6 + term4*term5*term3);
        }
    }
    return Tij;
}


double compute_Vab(STOGPrimitive gto1, STOGPrimitive gto2, double R[3]) {
    /*Calculate the nuclear attraction integral between two GTO primitives for a single nucleus*/
    double p = gto1.alpha + gto2.alpha;
    double P[3];
    P[0] = gto1.alpha * gto1.cords[0] / p + gto2.alpha * gto2.cords[0] / p;
    P[1] = gto1.alpha * gto1.cords[1] / p + gto2.alpha * gto2.cords[1] / p;
    P[2] = gto1.alpha * gto1.cords[2] / p + gto2.alpha * gto2.cords[2] / p;
    
    double val = 0.0;
    for (int t = 0; t <= gto1.nx + gto2.nx; t++) {
        for (int u = 0; u <= gto1.ny + gto2.ny; u++) {
            for (int v = 0; v <= gto1.nz + gto2.nz; v++) {
                double term1 = compute_E(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t);
                double term2 = compute_E(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u);
                double term3 = compute_E(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v);
                double term4 = compute_R(t, u, v, 0, p, P, R);
                val += term1 * term2 * term3 * term4;
            }
        }
    }
    return 2 * M_PI / p * val;
}

double compute_VijR(STOGOrbital orbital1, STOGOrbital orbital2, double R[3]) {        
    /*Compute the nuclear attraction integral between two STO-nG orbitals.*/
    int n = orbital1.n;
    int m = orbital2.n;
    double VijR = 0.0;
    
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOGPrimitive gto1 = orbital1.primitives[u];
            STOGPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.N * gto2.N * gto1.cc * gto2.cc;
            VijR += coeff * compute_Vab(gto1, gto2, R);
        }
    }
    return VijR;
}

double Vabcd(STOGPrimitive gto1, STOGPrimitive gto2, STOGPrimitive gto3, STOGPrimitive gto4) {
    /*Calculate the electron-electron repulsion integral between four GTO primitives. https://joshuagoings.com/2017/04/28/integrals/*/
    double p = gto1.alpha + gto2.alpha;
    double q = gto3.alpha + gto4.alpha;
    double alpha = p * q / (p + q);
    
    double P[3], Q[3];
    P[0] = gto1.alpha * gto1.cords[0] / p + gto2.alpha * gto2.cords[0] / p;
    P[1] = gto1.alpha * gto1.cords[1] / p + gto2.alpha * gto2.cords[1] / p;
    P[2] = gto1.alpha * gto1.cords[2] / p + gto2.alpha * gto2.cords[2] / p;
    Q[0] = gto3.alpha * gto3.cords[0] / q + gto4.alpha * gto4.cords[0] / q;
    Q[1] = gto3.alpha * gto3.cords[1] / q + gto4.alpha * gto4.cords[1] / q;
    Q[2] = gto3.alpha * gto3.cords[2] / q + gto4.alpha * gto4.cords[2] / q;
    
    double val = 0.0;
    for (int t = 0; t <= gto1.nx + gto2.nx; t++) {
        for (int u = 0; u <= gto1.ny + gto2.ny; u++) {
            for (int v = 0; v <= gto1.nz + gto2.nz; v++) {
                for (int tau = 0; tau <= gto3.nx + gto4.nx; tau++) {
                    for (int nu = 0; nu <= gto3.ny + gto4.ny; nu++) {
                        for (int phi = 0; phi <= gto3.nz + gto4.nz; phi++) {
                            double term1 = compute_E(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t);
                            double term2 = compute_E(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u);
                            double term3 = compute_E(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v);
                            double term4 = compute_E(gto3.cords[0], gto4.cords[0], gto3.alpha, gto4.alpha, gto3.nx, gto4.nx, tau);
                            double term5 = compute_E(gto3.cords[1], gto4.cords[1], gto3.alpha, gto4.alpha, gto3.ny, gto4.ny, nu);
                            double term6 = compute_E(gto3.cords[2], gto4.cords[2], gto3.alpha, gto4.alpha, gto3.nz, gto4.nz, phi);
                            double term7 = pow(-1, tau + nu + phi);
                            double term8 = compute_R(t + tau, u + nu, v + phi, 0, alpha, P, Q);
                            val += term1 * term2 * term3 * term4 * term5 * term6 * term7 * term8;
                        }
                    }
                }
            }
        }
    }
    double prefactor = 2 * pow(M_PI, 2.5) / (p * q * sqrt(p + q));
    return prefactor * val;
}

double compute_Vijkl(STOGOrbital orbital1, STOGOrbital orbital2, STOGOrbital orbital3, STOGOrbital orbital4) {
    /*Calculate the Coulomb repulsion integral between four STO-nG orbitals*/
    double Vijkl = 0.0;
    
    for (int m = 0; m < orbital1.n; m++) {
        for (int n = 0; n < orbital2.n; n++) {
            for (int u = 0; u < orbital3.n; u++) {
                for (int v = 0; v < orbital4.n; v++) {
                    STOGPrimitive gto1 = orbital1.primitives[m];
                    STOGPrimitive gto2 = orbital2.primitives[n];
                    STOGPrimitive gto3 = orbital3.primitives[u];
                    STOGPrimitive gto4 = orbital4.primitives[v];
                    
                    double norms = gto1.N * gto2.N * gto3.N * gto4.N;
                    double coefs = gto1.cc * gto2.cc * gto3.cc * gto4.cc;
                    Vijkl += norms * coefs * Vabcd(gto1, gto2, gto3, gto4);
                }
            }
        }
    }
    return Vijkl;
}
