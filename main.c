#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define _USE_MATH_DEFINES
#define N 256 // sim size
#define NUM_PARTICLES 20000 // number of particles

int i, j;

int mass_deposition_cic(double* position, double* target, int num_points) {

    for (int i = 0; i < num_points; ++i) { // iterate over all points in space

        double x = position[2 * i]; // x position of i-th particle
        double y = position[(2 * i) + 1]; // y position of i-th particle

        // find integer cell in which the particle lies
        int idx0 = (int)floor(x); // integer x cell of i-th particle
        int idx1 = (int)floor(y); // integer y cell of i-th particle

        double f0 = x - idx0; // fractional part in x
        double f1 = y - idx1; // fractional part in y

        // apply periodic boundary conditions
        int i0  = ((idx0) % N + N) % N;
        int i1  = ((idx1) % N + N) % N;
        int i0p = ((idx0 + 1) % N + N) % N;
        int i1p = ((idx1 + 1) % N + N) % N;

        // deposit onto mass density array
        target[i0 * N + i1] += (1 - f0) * (1 - f1);
        target[i0p *N + i1] += f0 * (1 - f1);
        target[i0 * N + i1p] += (1 - f0) * f1;
        target[i0p *N + i1p] += f0 * f1;
    }
};

// generate random number in range a to b
double uniform_random_double(double a, double b) {
    return a + (b - a) * (rand() / (double)RAND_MAX);
};

int main() {

    srand(time(NULL));

    // create and open file to store mass density data
    char const *file_name = "mass_density.dat";
    FILE *file = fopen(file_name, "w");
    
    // initialise density and potential arrays to zero
    double *density = calloc(N * N, sizeof(double));
    double *phi = calloc(N * N, sizeof(double));

    // initialise array to store particle positions
    double *particles = malloc(2 * NUM_PARTICLES * sizeof(double));

    // allocate random starting positions to the particles
    for (i = 0; i < NUM_PARTICLES; ++i) {
        particles[2 * i] = uniform_random_double(0.0, (double)N);
        particles[2 * i + 1] = uniform_random_double(0.0, (double)N);
    }

    // apply cic
    mass_deposition_cic(particles, density, NUM_PARTICLES);

    // print resulting mass distribution
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            fprintf(file, "%6.2f ", density[i * N + j]);
        };
        if (i != N - 1) fprintf(file, "\n");
    };

    fclose(file);

    // use command line to run gnuplot instructions from file
    system("gnuplot -persistent gnuplot_instructions.gp");

    return 0;
}