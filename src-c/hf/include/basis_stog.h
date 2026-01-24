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
} STOGPrimitive;

typedef struct {
    STOGPrimitive *primitives; // Array of GTO primitives
    int n; // Number of GTO primitives in the STO-nG
} STOGOrbital;

double compute_N(double alpha, int *n);

double compute_Sij(STOGOrbital orbital1, STOGOrbital orbital2); /*https://content.wolfram.com/sites/19/2012/02/Ho.pdf*/
double compute_Tij(STOGOrbital orbital1, STOGOrbital orbital2);
double compute_VijR(STOGOrbital orbital1, STOGOrbital orbital2, double R[3]);
double compute_Vijkl(STOGOrbital orbital1, STOGOrbital orbital2, STOGOrbital orbital3, STOGOrbital orbital4);
#endif
