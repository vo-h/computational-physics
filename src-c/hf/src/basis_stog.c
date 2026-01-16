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
double compute_nx(double A, double B, double alpha, double beta, int ai, int bi, double R, double t) {
    /*Compute the nuclear attraction integral between two STO-nG orbitals by multiplying the integrals along each axis and summing contributions from all GTO primitives.*/
    double P = (alpha * A + beta * B) / (alpha + beta);
    if (ai < 0 || bi < 0) {return 0.0;}
    if (ai == 0 && bi == 0) {return 1;}
    if (ai == 1 && bi == 0) {return -(A-P) - pow(t, 2)*(P-R);}
    
    if (bi == 0) {
        double term1 = -(A-P);
        double term2 = -pow(t, 2) * (P-R) * compute_nx(A, B, alpha, beta, ai-1, 0, R, t);
        double term3 = (ai-1) / (2 * (alpha + beta)) * (1-pow(t,2)) * compute_nx(A, B, alpha, beta, ai-2, 0, R, t);
        return term1 + term2 + term3;
    }

    double term1 = compute_nx(A, B, alpha, beta, ai+1, bi-1, R, t);
    double term2 = (A-B)*compute_nx(A, B, alpha, beta, ai, bi-1, R, t);
    return term1 + term2;
}

double integrand(STOPrimitive gto1, STOPrimitive gto2, double *R, double t) {
    t = 0.5*(t+1);
    double term1 = (gto1.alpha + gto2.alpha)*pow(t,2);
    double term2x = ((gto1.alpha * gto1.cords[0] + gto2.alpha * gto2.cords[0]) / (gto1.alpha + gto2.alpha)) - R[0];
    double term2y = ((gto1.alpha * gto1.cords[1] + gto2.alpha * gto2.cords[1]) / (gto1.alpha + gto2.alpha)) - R[1];
    double term2z = ((gto1.alpha * gto1.cords[2] + gto2.alpha * gto2.cords[2]) / (gto1.alpha + gto2.alpha)) - R[2];
    double term3 = term2x*term2x + term2y*term2y + term2z*term2z;
    double v_x = compute_nx(gto1.cords[0], gto2.cords[0], gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, R[0], t);
    double v_y = compute_nx(gto1.cords[1], gto2.cords[1], gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, R[1], t);
    double v_z = compute_nx(gto1.cords[2], gto2.cords[2], gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, R[2], t);
    return 0.5 * exp(-term1*term3) * v_x * v_y * v_z;
}

double abscissa(int n, int i) {
    /*Calculate the abscissa for the i-th point in an n-point Chebyshev quadrature.*/
    double term0 = i * M_PI / (n + 1);
    double term1 = (n+1-2*i) / (double)(n+1);
    double term2 = 1 + (2/3.0)*pow(sin(term0), 2);
    double term3 = cos(term0);
    double term4 = sin(term0);
    return term1 + (2/M_PI) * term2 * term3 * term4;
}

double omega(int n, int i) {
    /*Calculate the weight for the i-th point in an n-point Chebyshev quadrature.*/
    double term0 = i * (M_PI / (n + 1));
    return (16 / (3.0 * (n+1))) * pow(sin(term0), 4);
}

double int_chebyshev(double eps, int m, STOPrimitive gto1, STOPrimitive gto2, double *R) {
    /*Perform numerical integration of a function f using Chebyshev quadrature with m points and a specified tolerance eps.*/
    double err = 10.0;
    int n = 3;
    double c0 = cos(M_PI/6.0);
    double s0 = sin(M_PI/6.0);
    double c1 = s0;
    double s1 = c0;
    double q = (integrand(gto1, gto2, R, abscissa(2,1)) + integrand(gto1, gto2, R, -abscissa(2,1))) * omega(2,1);
    double p = integrand(gto1, gto2, R, 0.0);
    double chp = q+p;
    int j = 0;
    double c;
    double s;
    double xp;

    while ((err > eps) && (2*n*(1-j) + j*4*n/3.0 - 1 <= m)) {
        j = 1 - j;
        c1 = j * c1 + (1-j) * c0;
        s1 = j * s1 + (1-j) * s0;
        c0 = j * c0 + (1-j) * sqrt((1+c0) * 0.5);
        s0 = j * s0 + (1-j) * s0 / (c0+c0);
        c = c0;
        s = s0;

        for (int i = 1; i < n; i+=2) {
            xp = 1 + (2/(3*M_PI)) * s * c * (double)(3+2*s*s) - (i/ (double)(n));
            if (ceil(3.0*(i+j+j)/3.0) > i+j) {
                chp += (integrand(gto1, gto2, R, -xp) + integrand(gto1, gto2, R, xp)) * pow(s, 4);
            }
            xp = s;
            s = s*c1 + c*s1;
            c = c*c1 - xp*s1;
        }

        n = (1+j) * n;
        p = p + (1-j) * (chp-q);
        err = 16 * fabs((1-j)*(q-3*p/2.0) + j * (chp-2*q)) / (double)(3*n);
        q = (1-j) * q + j * chp;
    }
    return 16 * q / (double)(3*n);
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
            double coeff = gto1.cc * gto1.N * gto2.cc * gto2.N;
            double r2 = pow(gto1.cords[0] - gto2.cords[0], 2) + pow(gto1.cords[1] - gto2.cords[1], 2) + pow(gto1.cords[2] - gto2.cords[2], 2);
            double E_AB = exp(-gto1.alpha * gto2.alpha * r2 / (gto1.alpha + gto2.alpha));
            double prefactor = 2 * M_PI / (gto1.alpha + gto2.alpha);
            VijR += coeff * prefactor * E_AB * int_chebyshev(1e-10, 50000, gto1, gto2, R);
        }
    }
    return VijR;
}
