#include <math.h>
#include <stdio.h>

double cdiff1_1D(double (*func)(double), double x, double h) {
    /*Calculate the numerical gradient of a function at a given point using central difference.*/
    return (func(x + h) - func(x - h)) / (2 * h);
}

double cdiff1_2D(double (*func)(double, double), double x, double y, double h) {
    /*Calculate the numerical gradient of a function at a given point in 2D using central difference.*/
    double dfdx = (func(x + h, y) - func(x - h, y)) / (2 * h);
    double dfdy = (func(x, y + h) - func(x, y - h)) / (2 * h);
    return sqrt(dfdx*dfdx + dfdy*dfdy); // Return the magnitude of the gradient
}

double cdiff1_3D(double (*func)(double, double, double), double x, double y, double z, double h) {
    /*Calculate the numerical gradient of a function at a given point in 3D using central difference.*/
    double dfdx = (func(x + h, y, z) - func(x - h, y, z)) / (2 * h);
    double dfdy = (func(x, y + h, z) - func(x, y - h, z)) / (2 * h);
    double dfdz = (func(x, y, z + h) - func(x, y, z - h)) / (2 * h);
    return sqrt(dfdx*dfdx + dfdy*dfdy + dfdz*dfdz); // Return the magnitude of the gradient
}


double cdiff2_1D(double (*func)(double), double x, double h) {
    /*Calculate the numerical second derivative of a function at a given point using central difference.*/
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
}

double cdiff2_2D(double (*func)(double, double), double x, double y, double h) {
    /*Calculate the numerical second derivative of a function at a given point in 2D using central difference.*/
    double d2fdx2 = (func(x + h, y) - 2 * func(x, y) + func(x - h, y)) / (h * h);
    double d2fdy2 = (func(x, y + h) - 2 * func(x, y) + func(x, y - h)) / (h * h);
    return d2fdx2 + d2fdy2; // Return the Laplacian
}

double cdiff2_3D(double (*func)(double, double, double), double x, double y, double z, double h) {
    /*Calculate the numerical second derivative of a function at a given point in 3D using central difference.*/
    double d2fdx2 = (func(x + h, y, z) - 2 * func(x, y, z) + func(x - h, y, z)) / (h * h);
    double d2fdy2 = (func(x, y + h, z) - 2 * func(x, y, z) + func(x, y - h, z)) / (h * h);
    double d2fdz2 = (func(x, y, z + h) - 2 * func(x, y, z) + func(x, y, z - h)) / (h * h);
    return d2fdx2 + d2fdy2 + d2fdz2; // Return the Laplacian
}