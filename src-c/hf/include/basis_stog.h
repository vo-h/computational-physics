#ifndef BASIS_STOG_H // Use include guards to prevent multiple inclusions
#define BASIS_STOG_H

typedef struct {
    double cords[3]; // Coordinates of the center of the GTO
    double alpha; // Exponents for the GTO primitives
    double cc; // Coefficients for the linear combination of GTOs
    double N; // Normalization constant for the GTO primitive
    int nx; // Angular momentum quantum numbers along x-axis
    int ny; // Angular momentum quantum numbers along y-axis
    int nz; // Angular momentum quantum numbers along z-axis
} STOPrimitive;

typedef struct {
    STOPrimitive *primitives; // Array of GTO primitives
    int n; // Number of GTO primitives in the STO-nG
} STOOrbital ;

double compute_N(double alpha, int *n);

double compute_Sij(STOOrbital orbital1, STOOrbital orbital2); /*https://content.wolfram.com/sites/19/2012/02/Ho.pdf*/
double compute_Tij(STOOrbital orbital1, STOOrbital orbital2);
double compute_Vij(STOOrbital orbital1, STOOrbital orbital2);

double compute_Jijkl(STOOrbital orbital1, STOOrbital orbital2, STOOrbital orbital3, STOOrbital orbital4);
double compute_Kijkl(STOOrbital orbital1, STOOrbital orbital2, STOOrbital orbital3, STOOrbital orbital4);

#endif
