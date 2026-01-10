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

double interate_ss(double *COORS1, double *COORS2, double alpha1, double alpha2) {
    /*Calculate the overlap integral between two s-type GTO primitives.*/
    double pi = M_PI;
    double prefactor = sqrt(M_PI / (alpha1 + alpha2));
    double a = alpha1 + alpha2;
    
    double bx = 2*alpha1*COORS1[0] + 2*alpha2*COORS2[0];
    double by = 2*alpha1*COORS1[1] + 2*alpha2*COORS2[1];
    double bz = 2*alpha1*COORS1[2] + 2*alpha2*COORS2[2];

    double cx = alpha1 * pow(COORS1[0], 2) + alpha2 * pow(COORS2[0], 2);
    double cy = alpha1 * pow(COORS1[1], 2) + alpha2 * pow(COORS2[1], 2);
    double cz = alpha1 * pow(COORS1[2], 2) + alpha2 * pow(COORS2[2], 2);

    double resx = exp((bx*bx / (4*a)) - cx);
    double resy = exp((by*by / (4*a)) - cy);
    double resz = exp((bz*bz / (4*a)) - cz);
    return pow(prefactor, 3) * resx * resy * resz;
}

double integrate_sp(double *COORS_s, double *COORS_p, int *n_s, int *n_p, double alpha_s, double alpha_p) {
    /*Calculate the overlap integral between an s-type and a p-type GTO primitive.
      Formula: (alpha_s * (X_p - X_s)) / (alpha_s + alpha_p) * integral_ss
    */
    
    double center_s;
    double center_p;
    if (n_p[0] == 1) {
        center_s = COORS_s[0];
        center_p = COORS_p[0];
    } else if (n_p[1] == 1) {
        center_s = COORS_s[1];
        center_p = COORS_p[1];
    } else if (n_p[2] == 1) {
        center_s = COORS_s[2];
        center_p = COORS_p[2];
    } else {
        return 0.0; // Not a valid s-p combination
    }

    double prefactor = (alpha_s * (center_p - center_s)) / (alpha_s + alpha_p);
    double overlap_ss = interate_ss(COORS_s, COORS_p, alpha_s, alpha_p);
    return -1 * prefactor * overlap_ss;
}

double integrate_pp(double *COORS1, double *COORS2, int *n1, int *n2, double alpha1, double alpha2) {
    /*Calculate the overlap integral between two p-type GTO primitives.
      
      For same orbital type (e.g., px-px):
          [1/(2(alpha1+alpha2)) - alpha1*alpha2*(Xi-Xj)^2/(alpha1+alpha2)^2] * integral_ss
      
      For orthogonal orbitals (e.g., px-py):
          [alpha1*(Xi-Xj)/(alpha1+alpha2)] * [alpha1*(Yi-Yj)/(alpha1+alpha2)] * integral_ss
    */
    
    double alpha_sum = alpha1 + alpha2;
    double ss_integral = interate_ss(COORS1, COORS2, alpha1, alpha2);
    
    // Check if same orbital type (e.g., px-px, py-py, pz-pz)
    if ((n1[0] == 1 && n2[0] == 1) || (n1[1] == 1 && n2[1] == 1) || (n1[2] == 1 && n2[2] == 1)) {
        // Get the relevant coordinate difference
        double diff;
        if (n1[0] == 1) {
            diff = COORS1[0] - COORS2[0];
        } else if (n1[1] == 1) {
            diff = COORS1[1] - COORS2[1];
        } else { // n1[2] == 1
            diff = COORS1[2] - COORS2[2];
        }
        
        double prefactor = (1.0 / (2.0 * alpha_sum)) - (alpha1 * alpha2 * diff * diff) / (alpha_sum * alpha_sum);
        return prefactor * ss_integral;
    }
    
    // Orthogonal orbitals (e.g., px-py, px-pz, py-pz)
    else {
        double factor = 1.0;
        
        // Check each dimension and accumulate factors
        if (n1[0] == 1) {
            double diff_x = COORS1[0] - COORS2[0];
            factor *= (alpha1 * diff_x) / alpha_sum;
        }
        if (n1[1] == 1) {
            double diff_y = COORS1[1] - COORS2[1];
            factor *= (alpha1 * diff_y) / alpha_sum;
        }
        if (n1[2] == 1) {
            double diff_z = COORS1[2] - COORS2[2];
            factor *= (alpha1 * diff_z) / alpha_sum;
        }
        
        if (n2[0] == 1 && n1[0] == 0) {
            double diff_x = COORS1[0] - COORS2[0];
            factor *= (alpha2 * diff_x) / alpha_sum;
        }
        if (n2[1] == 1 && n1[1] == 0) {
            double diff_y = COORS1[1] - COORS2[1];
            factor *= (alpha2 * diff_y) / alpha_sum;
        }
        if (n2[2] == 1 && n1[2] == 0) {
            double diff_z = COORS1[2] - COORS2[2];
            factor *= (alpha2 * diff_z) / alpha_sum;
        }
        
        return factor * ss_integral;
    }
}

// double Sij(double *COORS1, double *COORS2, double *cc1, double *cc2, double *alpha1, double *alpha2, int *n1, int *n2, int n) {
//     /*Calculate the overlap integral S_ij between two GTOs centered at (X1, Y1, Z1) and (X2, Y2, Z2) with exponents alpha1 and alpha2 and angular momentum quantum numbers (nx1, ny1, nz1) and (nx2, ny2, nz2).*/
//     // This function is not implemented in this code snippet but can be implemented using the Gaussian product theorem and analytical formulas for overlap integrals.
//     return 0.0; // Placeholder return value
// }

