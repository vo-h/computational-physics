#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Chebyshev quadrature for numerical integration. */
double abscissa(int n, int i) {
    /*Calculate the abscissa for the i-th point in an n-point Chebyshev quadrature.*/
    double term0 = i * M_PI / (n + 1);
    double term1 = (n+1-2*i) / (n+1);
    double term2 = 1 + (2/3)*pow(sin(term0), 2);
    double term3 = cos(term0);
    double term4 = sin(term0);
    return term1 + (2/M_PI) * term2 * term3 * term4;
}

double omega(int n, int i) {
    /*Calculate the weight for the i-th point in an n-point Chebyshev quadrature.*/
    double term0 = i + M_PI / (n + 1);
    return (16/3*(n+1)) * pow(sin(term0), 4);
}

double int_chebyshev(double eps, int m, double (*f)(double)) {
    /*Perform numerical integration of a function f using Chebyshev quadrature with m points and a specified tolerance eps.*/
    int err = 10;
    int n = 3;
    double c0 = cos(M_PI/6);
    double s0 = sin(M_PI/2);
    double c1 = s0;
    double s1 = c0;
    double q = f(abscissa(2,1)) + f(-abscissa(2,1)) * omega(2,1);
    double p = f(0.0);
    double chp = q+p;
    int j = 0;
    double c;
    double s;
    double xp;

    while ((err > eps) && (2*n*(1-j) + j*4*n/3 - 1 < m)) {
        j = 1 - j;
        c1 = j * c1 + (1-j) * c0;
        s1 = j * s1 + (1-j) * s0;
        c0 = j * c0 + (1-j) * sqrt(0.5*(1+c0));
        s0 = j * s0 + (1-j) * s0 / (c0+c0);
        c = c0;
        s = s0;

        for (int i = 1; i <= n; i+=2) {
            xp = 1 + (2/(3*M_PI)) * s * c * (3+2*s*s) - i/n;
            if (ceil(3*(i+j+j)/3) > i+j) {
                chp += (f(-xp) + f(xp)) * pow(s, 4);
            }
            xp = s;
            s = s*c1 + c*s1;
            c = c*c1 - xp*s1;
        }

        n = (1+j) * n;
        p = p + (1-j) * (chp-q);
        err = 16 * fabs((1-j)*(q-3*p/2) + j * (chp-2*q)) / (3*n);
        q = (1-j) * q + j * chp;
    }
    return 16 * q / (3*n);
}