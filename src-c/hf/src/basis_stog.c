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

/*Function to compute the overlap integral matrix: https://content.wolfram.com/sites/19/2012/02/Ho.pdf*/
double compute_sx(double A, double B, double alpha, double beta, int ai, int bi) {
    /*Compute the integral of the product of two GTOs along a single axis (x, y, or z) using recursion.*/

    double P = (alpha * A + beta * B) / (alpha + beta);

    /*Base cases for recursion*/
    if (ai < 0 || bi < 0) {return 0.0;}
    if (ai == 0 && bi == 0) {return 1;}
    if (ai == 1 && bi == 0) {return -(A-P);}

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
            double contribution = coeff * prefactor * E_AB * sx * sy * sz;
            Sij += contribution;
        }
    }
    return Sij;
}

/*https://content.wolfram.com/sites/19/2013/01/Ho_Kinetic.pdf*/
double compute_tx(double A, double B, double alpha, double beta, int ai, int bi) {
    /*Compute the kinetic energy integral between two GTOs along a single axis using recursion.*/

    double P = (alpha * A + beta * B) / (alpha + beta);

    /*Base cases for recursion*/
    if (ai < 0 || bi < 0) {return 0.0;}
    if (ai == 0 && bi == 0) {return 2*alpha*beta*compute_sx(A, B, alpha, beta, 1, 1);}
    if (bi == 0) {return -ai*beta*compute_sx(A, B, alpha, beta, ai-1, 1) + 2*alpha*beta*compute_sx(A, B, alpha, beta, ai+1, 1);}
    if (ai == 0) {return -bi*alpha*compute_sx(A, B, alpha, beta, 1, bi-1) + 2*alpha*beta*compute_sx(A, B, alpha, beta, 1, bi+1);}

    /*Index recursion: eq. 11*/
    double term1 = ai*bi*compute_sx(A, B, alpha, beta, ai-1, bi-1);
    double term2 = 2*ai*beta*compute_sx(A, B, alpha, beta, ai-1, bi+1);
    double term3 = 2*bi*alpha*compute_sx(A, B, alpha, beta, ai+1, bi-1);
    double term4 = 4*alpha*beta*compute_sx(A, B, alpha, beta, ai+1, bi+1);
    return 0.5*(term1 - term2 - term3 + term4);
}

double compute_Tij(STOOrbital orbital1, STOOrbital orbital2) {
    /*Compute the kinetic energy integral between two STO-nG orbitals by multiplying the integrals along each axis and summing contributions from all GTO primitives.*/

    int n = orbital1.n; 
    int m = orbital2.n;
    double Tij = 0.0;
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOPrimitive gto1 = orbital1.primitives[u];
            STOPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.cc * gto1.N * gto2.cc * gto2.N;
            double r2 = pow(gto1.cords[0] - gto2.cords[0], 2) + pow(gto1.cords[1] - gto2.cords[1], 2) + pow(gto1.cords[2] - gto2.cords[2], 2);
            double E_AB = exp(-gto1.alpha * gto2.alpha * r2 / (gto1.alpha + gto2.alpha));
            double prefactor = pow(M_PI / (gto1.alpha + gto2.alpha), 1.5);

            double tx = compute_tx(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx);
            double ty = compute_tx(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny);
            double tz = compute_tx(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz);

            double sx = compute_sx(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx);
            double sy = compute_sx(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny);
            double sz = compute_sx(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz);
            Tij += coeff * prefactor * E_AB * (tx*sy*sz + sx*ty*sz + sx*sy*tz);
        }
    }
    return Tij;
}

/*https://content.wolfram.com/sites/19/2014/12/Ho_Nuclear.pdf*/
double compute_boys(int n, double T) {
    /*Compute the Boys function using the confluent hypergeometric function 1F1, which is used in the evaluation of nuclear attraction integrals.*/
    return gsl_sf_hyperg_1F1(n+0.5, n+1.5, -T) / (2.0*n + 1.0);
}
double compute_E(int ai, int bi, int t, double A, double B, double alpha, double beta) {
    /*Compute the nuclear attraction integral between two GTOs along a single axis using recursion.*/
    double p = alpha + beta;
    double q = alpha * beta / p;
    double Qx = A - B;

    if (ai < 0 || t > ai+bi) {return 0.0;}
    if (ai == 0 && bi == 0 && t == 0) {return exp(-q * Qx * Qx);}
    if (bi == 0) {
        double term1 = (1/(2*p)) * compute_E(ai-1, bi, t-1, A, B, alpha, beta);
        double term2 = (q * Qx / alpha) * compute_E(ai-1, bi, t, A, B, alpha, beta);
        double term3 = (t+1) * compute_E(ai-1, bi, t+1, A, B, alpha, beta);
        return term1 - term2 + term3;
    }

    double term1 = (1/(2*p)) * compute_E(ai, bi-1, t-1, A, B, alpha, beta);
    double term2 = (q * Qx / beta) * compute_E(ai, bi-1, t, A, B, alpha, beta);
    double term3 = (t+1) * compute_E(ai, bi-1, t+1, A, B, alpha, beta);
    return term1 + term2 + term3;
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
    double T = p * RPC * RPC;
    double val = 0.0;
    
    if (t == 0 && u == 0 && v == 0) {
        val += pow(-2*p, n) * compute_boys(n, T);
    } else if (t == 0 && u == 0) {
        if (v > 1) {
            val += (v-1) * compute_R(t, u, v-2, n+1, p, P, C);
        }
        val += (P[2]-C[2]) * compute_R(t, u, v-1, n+1, p, P, C);
    } else if (t == 0) {
        if (u > 1) {
            val += (u-1) * compute_R(t, u-2, v, n+1, p, P, C);
        }
        val += (P[1]-C[1]) * compute_R(t, u-1, v, n+1, p, P, C);
    } else {
        if (t > 1) {
            val += (t-1) * compute_R(t-2, u, v, n+1, p, P, C);
        }
        val += (P[0]-C[0]) * compute_R(t-1, u, v, n+1, p, P, C);
    }
    return val;
}

double compute_Pi(STOPrimitive gto1, STOPrimitive gto2, int component) {
    /*Calculate the weighted center for the integral in a given component (0=x, 1=y, 2=z)*/
    double alpha_sum = gto1.alpha + gto2.alpha;
    return (gto1.alpha * gto1.cords[component] + gto2.alpha * gto2.cords[component]) / alpha_sum;
}

double compute_nuclear_attraction(STOPrimitive gto1, STOPrimitive gto2, double R[3]) {
    /*Calculate the nuclear attraction integral between two GTO primitives for a single nucleus*/
    double p = gto1.alpha + gto2.alpha;
    double P[3];
    P[0] = compute_Pi(gto1, gto2, 0);
    P[1] = compute_Pi(gto1, gto2, 1);
    P[2] = compute_Pi(gto1, gto2, 2);
    
    double val = 0.0;
    for (int t = 0; t <= gto1.nx + gto2.nx; t++) {
        for (int u = 0; u <= gto1.ny + gto2.ny; u++) {
            for (int v = 0; v <= gto1.nz + gto2.nz; v++) {
                double term1 = compute_E(gto1.nx, gto2.nx, t, gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha);
                double term2 = compute_E(gto1.ny, gto2.ny, u, gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha);
                double term3 = compute_E(gto1.nz, gto2.nz, v, gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha);
                double term4 = compute_R(t, u, v, 0, p, P, R);
                val += term1 * term2 * term3 * term4;
            }
        }
    }
    return 2 * M_PI / p * val;
}

double compute_VijR(STOOrbital orbital1, STOOrbital orbital2, double R[3]) {        
    /*Compute the nuclear attraction integral between two STO-nG orbitals by multiplying the integrals along each axis and summing contributions from all GTO primitives.*/
    int n = orbital1.n;
    int m = orbital2.n;
    double VijR = 0.0;
    
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < m; v++) {
            STOPrimitive gto1 = orbital1.primitives[u];
            STOPrimitive gto2 = orbital2.primitives[v];
            double coeff = gto1.N * gto2.N * gto1.cc * gto2.cc;
            VijR += coeff * compute_nuclear_attraction(gto1, gto2, R);
        }
    }
    return VijR;
}

double compute_repulsion(STOPrimitive gto1, STOPrimitive gto2, STOPrimitive gto3, STOPrimitive gto4) {
    /*Calculate the electron-electron repulsion integral between four GTO primitives. https://joshuagoings.com/2017/04/28/integrals/*/
    double p = gto1.alpha + gto2.alpha;
    double q = gto3.alpha + gto4.alpha;
    double alpha = p * q / (p + q);
    
    double P[3], Q[3];
    P[0] = compute_Pi(gto1, gto2, 0);
    P[1] = compute_Pi(gto1, gto2, 1);
    P[2] = compute_Pi(gto1, gto2, 2);
    Q[0] = compute_Pi(gto3, gto4, 0);
    Q[1] = compute_Pi(gto3, gto4, 1);
    Q[2] = compute_Pi(gto3, gto4, 2);
    
    double val = 0.0;
    for (int t = 0; t <= gto1.nx + gto2.nx; t++) {
        for (int u = 0; u <= gto1.ny + gto2.ny; u++) {
            for (int v = 0; v <= gto1.nz + gto2.nz; v++) {
                for (int tau = 0; tau <= gto3.nx + gto4.nx; tau++) {
                    for (int nu = 0; nu <= gto3.ny + gto4.ny; nu++) {
                        for (int phi = 0; phi <= gto3.nz + gto4.nz; phi++) {
                            double term1 = compute_E(gto1.nx, gto2.nx, t, gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha);
                            double term2 = compute_E(gto1.ny, gto2.ny, u, gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha);
                            double term3 = compute_E(gto1.nz, gto2.nz, v, gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha);
                            double term4 = compute_E(gto3.nx, gto4.nx, tau, gto3.cords[0], gto4.cords[0], gto3.alpha, gto4.alpha);
                            double term5 = compute_E(gto3.ny, gto4.ny, nu, gto3.cords[1], gto4.cords[1], gto3.alpha, gto4.alpha);
                            double term6 = compute_E(gto3.nz, gto4.nz, phi, gto3.cords[2], gto4.cords[2], gto3.alpha, gto4.alpha);
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

double compute_Vijkl(STOOrbital orbital1, STOOrbital orbital2, STOOrbital orbital3, STOOrbital orbital4) {
    /*Calculate the Coulomb repulsion integral between four STO-nG orbitals*/
    double Vijkl = 0.0;
    
    for (int m = 0; m < orbital1.n; m++) {
        for (int n = 0; n < orbital2.n; n++) {
            for (int u = 0; u < orbital3.n; u++) {
                for (int v = 0; v < orbital4.n; v++) {
                    STOPrimitive gto1 = orbital1.primitives[m];
                    STOPrimitive gto2 = orbital2.primitives[n];
                    STOPrimitive gto3 = orbital3.primitives[u];
                    STOPrimitive gto4 = orbital4.primitives[v];
                    
                    double norms = gto1.N * gto2.N * gto3.N * gto4.N;
                    double coefs = gto1.cc * gto2.cc * gto3.cc * gto4.cc;
                    Vijkl += norms * coefs * compute_repulsion(gto1, gto2, gto3, gto4);
                }
            }
        }
    }
    return Vijkl;
}
