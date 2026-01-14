#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include "basis_stog.h"
#include <unistd.h>

int get_atomic_number(char *symbol) {
    if (strcmp(symbol, "H") == 0) {return 1;}
    if (strcmp(symbol, "He") == 0) {return 2;}
    if (strcmp(symbol, "Li") == 0) {return 3;}
    if (strcmp(symbol, "Be") == 0) {return 4;}
    if (strcmp(symbol, "B") == 0) {return 5;}
    if (strcmp(symbol, "C") == 0) {return 6;}
    if (strcmp(symbol, "N") == 0) {return 7;}
    if (strcmp(symbol, "O") == 0) {return 8;}
    if (strcmp(symbol, "F") == 0) {return 9;}
    if (strcmp(symbol, "Ne") == 0) {return 10;}
    if (strcmp(symbol, "Na") == 0) {return 11;}
    if (strcmp(symbol, "Mg") == 0) {return 12;}
    if (strcmp(symbol, "Al") == 0) {return 13;}
    if (strcmp(symbol, "Si") == 0) {return 14;}
    if (strcmp(symbol, "P") == 0) {return 15;}
    if (strcmp(symbol, "S") == 0) {return 16;}
    if (strcmp(symbol, "Cl") == 0) {return 17;}

}
/*https://www.basissetexchange.org/basis/sto-3g/format/gaussian94/?version=1&elements=8&uncontract_spdf=true*/

typedef struct  {
    char atom[3]; // Element symbol (e.g., "H", "C", "O")
    int Z; // Atomic number
    double coords[3]; // Atomic coordinates (x, y, z)
    STOOrbital *orbitals; // Array of STO-nG orbitals for the atom
    int num_orbitals; // Number of orbitals for the atom
} Atom;

Atom parse_atom(char *atom, int Z, double x, double y, double z, int num_orbitals) {
    Atom temp = {
        .atom = *atom,
        .Z = Z,
        .coords = {x, y, z},
        .orbitals = NULL,
        .num_orbitals = num_orbitals,
    };
    return temp;
}