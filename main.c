#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define _USE_MATH_DEFINES
#define N 256 // sim size
#define NUM_PARTICLES 100000 // number of particles

typedef struct{

    double** positions;
    int num_particles;

} Particles2d;

Particles2d* create_particles(int num_particles) {

    Particles2d* p = malloc(sizeof(Particles2d));
    p -> num_particles = num_particles;

    p -> positions = malloc(num_particles * sizeof(double*));
    for (int i = 0; i < num_particles; ++i) {
        p -> positions[i] = malloc(2 * sizeof(double));
    }

    return p;

}

void free_particles(Particles2d* p) {

    for (int i = 0; i < p -> num_particles; ++i) {
        free(p -> positions[i]);
    }

    free(p -> positions);
    free(p);

}

typedef struct {

    int grid_size;
    double* density_grid;

} CICDepositor;

CICDepositor* create_cic_depositor() {

    CICDepositor* cic = malloc(sizeof(CICDepositor));
    cic -> grid_size = N;
    cic -> density_grid = calloc(N * N, sizeof(double));
    return cic;

}

void free_cic_depositor(CICDepositor* cic) {

    free(cic -> density_grid);
    free(cic);

}

void deposit_mass(CICDepositor* cic, Particles2d* p) {

    for (int i = 0; i < p -> num_particles; ++i) {

        double x = p -> positions[i][0];
        double y = p -> positions[i][1];

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
        cic -> density_grid[i0 * N + i1] += (1 - f0) * (1 - f1);
        cic -> density_grid[i0p *N + i1] += f0 * (1 - f1);
        cic -> density_grid[i0 * N + i1p] += (1 - f0) * f1;
        cic -> density_grid[i0p *N + i1p] += f0 * f1;
    }
}

typedef struct {

    double** data;
    int rows, cols;

} Interpolation2d;

Interpolation2d* create_Interpolation(double** data, int rows, int cols) {

    Interpolation2d* interp = (Interpolation2d*)malloc(sizeof(Interpolation2d));

    interp -> data = data;
    interp -> rows = rows;
    interp -> cols = cols;

    return interp;
}

double* interp2d_call(Interpolation2d* interp, Particles2d* p) {

    double* result = (double*)malloc((p -> num_particles) * sizeof(double));

    for (int idx = 0; idx < p -> num_particles; ++idx) {

        double x = p -> positions[idx][0];
        double y = p -> positions[idx][1];

        int i0 = ((int)floor(x)) % (interp -> rows);
        int j0 = ((int)floor(x)) % (interp -> cols);
        int i1 = ((int)ceil(y)) % (interp -> rows);
        int j1 = ((int)ceil(y)) % (interp -> cols);

        double xm0 = x - floor(x);
        double xm1 = y - floor(y);
        double xn0 = 1.0 - xm0;
        double xn1 = 1.0 - xm1;

        double f1 = interp -> data[i0][j0];
        double f2 = interp -> data[i1][j0];
        double f3 = interp -> data[i0][j1];
        double f4 = interp -> data[i1][j1];

        result[idx] = (f1 * xn0 * xn1) + (f2 * xm0 * xn1) + (f3 * xn0 * xm1) + (f4 * xm0 * xm1);
    }

    return result;

}

double** second_order_gradient(double** F, int rows, int cols, int axis) {

    double** G = (double**)malloc(rows * sizeof(double*));
    
    for (int i =  0; i < rows; ++i) {
        G[i] = (double*)malloc(cols * sizeof(double));
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {

            int i_m2 = (i - 2 + rows) % rows;
            int i_m1 = (i - 1 + rows) % rows;
            int i_p1 = (i + 1) % rows;
            int i_p2 = (i + 2) % rows;

            int j_m2 = (j - 2 + cols) % cols;
            int j_m1 = (j - 1 + cols) % cols;
            int j_p1 = (j + 1) % cols;
            int j_p2 = (j + 2) % cols;

            if (axis == 0) {
                G[i][j] = ((1.0 / 12) * F[i_m2][j]) - ((2.0 / 3) * F[i_m1][j]) + ((2.0 / 3) * F[i_p1][j]) - ((1.0 / 12) * F[i_p2][j]);
            }

            else {
                G[i][j] = ((1.0/ 12) * F[i][j_m2]) - ((2.0 / 3) * F[i][j_m1]) + ((2.0 / 3) * F[i][j_p1]) - ((1.0 / 12) * F[i][j_p2]);
            }
        }
    }

    return G;

}

// generate random number in range a to b
double uniform_random_double(double a, double b) {

    return a + (b - a) * (rand() / (double)RAND_MAX);

}

int main() {

    srand(time(NULL));

    // create and open file to store mass density data
    char const *file_name = "mass_density.dat";
    FILE *file = fopen(file_name, "w");

    // create particles
    Particles2d* particles = create_particles(NUM_PARTICLES);

    // assign random positions
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        particles -> positions[i][0] = uniform_random_double(0.0, (double)N);
        particles -> positions[i][1] = uniform_random_double(0.0, (double)N);
    }

    // create cic
    CICDepositor* cic = create_cic_depositor();

    // deposit mass onto density grid
    deposit_mass(cic, particles);

    // print resulting mass distribution
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(file, "%6.2f ", cic -> density_grid[i * N + j]);
        };
        if (i != N - 1) fprintf(file, "\n");
    };

    fclose(file);

    // check if gnuplot is installed 
    // if so use command line to run gnuplot instructions from file
    int status = system("gnuplot --version >null 2>&1");
    if (status == 0) {
        system("gnuplot -persistent gnuplot_instructions.gp");
    }
    else  {
        (printf("gnuplot is not installed \n"));
        return 1;
    }

    // cleanup
    free_particles(particles);
    free_cic_depositor(cic);
    
    return 0;
}