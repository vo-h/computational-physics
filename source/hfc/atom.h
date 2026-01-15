#pragma once
#define _USE_MATH_DEFINES // Required for M_PI on some compilers
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "basis_stog.h"
#include <unistd.h>
#include <curl/curl.h>

typedef struct  {
    char atom[3]; // Element symbol (e.g., "H", "C", "O")
    int Z; // Atomic number
    double coords[3]; // Atomic coordinates (x, y, z)
    STOOrbital *orbitals; // Array of STO-nG orbitals for the atom
    int num_orbitals; // Number of orbitals for the atom
    int num_gtos; // Number of GTO primitives per orbital
} Atom;

void replace_char(char* str, char find, char replace) {
    int i = 0;
    while (str[i] != '\0') { // Loop until the null terminator is reached
        if (str[i] == find) {
            str[i] = replace; // Replace the character if it matches
        }
        i++;
    }
}

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
    return -1;
}

// Structure to hold response data from cURL
typedef struct {
    char *data;
    size_t size;
} CurlResponse;

// Callback function for cURL to write data
static size_t write_callback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t total_size = size * nmemb;
    CurlResponse *response = (CurlResponse *)userp;
    
    char *ptr = (char *)realloc(response->data, response->size + total_size + 1);
    if (ptr == NULL) {
        fprintf(stderr, "Not enough memory for cURL response\n");
        return 0;
    }
    
    response->data = ptr;
    memcpy(&(response->data[response->size]), contents, total_size);
    response->size += total_size;
    response->data[response->size] = 0; // Null terminate
    
    return total_size;
}

// Function to fetch orbital data from URL and parse into STOOrbital structures
STOOrbital* fetch_orbital_data(const char *url, const char *atom_symbol, double coords[3], int *num_orbitals, int num_gtos) {
    CURL *curl;
    CURLcode res;
    CurlResponse response = {.data = NULL, .size = 0};
    
    // Initialize cURL
    curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "Failed to initialize cURL\n");
        *num_orbitals = 0;
        return NULL;
    }
    
    // Set cURL options
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&response);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "libcurl-agent/1.0");
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    
    // Perform the request
    res = curl_easy_perform(curl);
    
    if (res != CURLE_OK) {
        fprintf(stderr, "cURL failed: %s\n", curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        free(response.data);
        *num_orbitals = 0;
        return NULL;
    }
    
    curl_easy_cleanup(curl);
    
    if (response.data == NULL) {
        *num_orbitals = 0;
        return NULL;
    }
    
    // Debug: print raw response
    // printf("=== Basis set data received ===\n%s\n=== End of data ===\n", response.data);
    
    // First pass: count shells (which become orbitals)
    // S shell = 1 orbital, P shell = 3 orbitals (px, py, pz), D shell = 5 orbitals
    int orbital_count = 0;
    char *line = response.data;
    char *next_line;
    char *data_copy = strdup(response.data);
    
    line = data_copy;
    while ((next_line = strchr(line, '\n')) != NULL) {
        *next_line = '\0';
        if (strlen(line) > 0 && line[1] == ' ') {
            if (line[0] == 'S') {
                orbital_count += 1;
            }
            else if (line[0] == 'P') {
                orbital_count += 3;
            }
            else if (line[0] == 'D') {
                orbital_count += 5;
            }
        }
        line = next_line + 1;
    }
    free(data_copy);
    
    if (orbital_count == 0) {
        free(response.data);
        *num_orbitals = 0;
        return NULL;
    }
    
    // Allocate array for orbitals
    STOOrbital *orbitals = (STOOrbital *)malloc(orbital_count * sizeof(STOOrbital));
    if (orbitals == NULL) {
        fprintf(stderr, "Failed to allocate memory for orbitals\n");
        free(response.data);
        *num_orbitals = 0;
        return NULL;
    }
    
    // Second pass: parse shells and create orbitals

    line = response.data;
    char shell_type = '\0';
    int orbital_idx = 0;
    double alpha;
    double cc;
    while (1) {
        next_line = strchr(line, '\n');
        if (next_line) *next_line = '\0';
        
        // Detect shell type line (e.g., "S   3   1.00")
        if (strlen(line) > 0 &&  strchr("SPD", line[0]) && line[1] == ' ') {

            shell_type = line[0];
            if (shell_type == 'S') {
                STOPrimitive *primitives = (STOPrimitive *)malloc(num_gtos * sizeof(STOPrimitive));

                for (int i = 0; i < num_gtos; i++) {
                    // Move line pointer to next primitive line
                    line = next_line + 1;
                    next_line = strchr(line, '\n');
                    if (next_line) *next_line = '\0';

                    // Parse primitive line
                    replace_char(line, 'D', 'E'); // Replace Fortran-style exponent if present
                    sscanf(line, "%lf %lf", &alpha, &cc);
                    int n[3] = {0, 0, 0};
                    STOPrimitive primitive = {
                        .alpha = alpha,
                        .cc = cc,
                        .cords = {coords[0], coords[1], coords[2]},
                        .nx = 0,
                        .ny = 0,
                        .nz = 0,
                        .N = compute_N(alpha, n)
                    };
                    primitives[i] = primitive;
                }
                orbitals[orbital_idx].primitives = primitives;
                orbitals[orbital_idx].n = num_gtos;
                orbital_idx++;
            } else if (shell_type == 'P') {
                double *alphas = (double *)malloc(num_gtos * sizeof(double));
                double *ccs = (double *)malloc(num_gtos * sizeof(double));
                for (int i = 0; i < num_gtos; i++) {
                    // Move line pointer to next primitive line
                    line = next_line + 1;
                    next_line = strchr(line, '\n');
                    if (next_line) *next_line = '\0';
                    replace_char(line, 'D', 'E'); // Replace Fortran-style exponent if present
                    sscanf(line, "%lf %lf", &alpha, &cc);
                    alphas[i] = alpha;
                    ccs[i] = cc;
                }

                // Create 3 orbitals (px, py, pz) with same primitives but different angular momentum
                for (int p_idx = 0; p_idx < 3; p_idx++) {
                    STOPrimitive *primitives = (STOPrimitive *)malloc(num_gtos * sizeof(STOPrimitive));
                    for (int i = 0; i < num_gtos; i++) {
                        int n[3] = {0, 0, 0};
                        if (p_idx == 0) n[0] = 1;
                        else if (p_idx == 1) n[1] = 1;
                        else if (p_idx == 2) n[2] = 1;
                        STOPrimitive primitive = {
                            .alpha = alphas[i],
                            .cc = ccs[i],
                            .cords = {coords[0], coords[1], coords[2]},
                            .nx = n[0],
                            .ny = n[1],
                            .nz = n[2],
                            .N = compute_N(alphas[i], n)
                        };
                        primitives[i] = primitive;
                    }
                    orbitals[orbital_idx].primitives = primitives;
                    orbitals[orbital_idx].n = num_gtos;
                    orbital_idx++;
                }
            };
            line = next_line + 1;
            
        }
        if (next_line == NULL) break;
        line = next_line + 1;
        if (line[0] == '*') break;
    }
    
    free(response.data);
    *num_orbitals = orbital_idx;
    return orbitals;
}

Atom parse_atom(char *symbol, double x, double y, double z, int num_gtos) {

    char url[128];
    int Z = get_atomic_number(symbol);
    int res = snprintf(url, sizeof(url), "https://www.basissetexchange.org/basis/sto-%dg/format/gaussian94/?version=1&elements=%d&uncontract_spdf=true", num_gtos, Z);
    if (res >= sizeof(url)) {
        printf("Warning: Output was truncated.\n");
    }

    Atom atom = {
        .atom = *symbol,
        .Z = Z,
        .coords = {x, y, z},
        .num_gtos = num_gtos,
    };
    atom.orbitals = fetch_orbital_data(url, symbol, atom.coords, &atom.num_orbitals, num_gtos);
    
    return atom;
}