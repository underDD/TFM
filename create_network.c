#include "head.h"

int main(){

    // ? ############### VARIABLES ###############

    double l_mean, l_sigma; 
    double d_mean, d_sigma;

    int N_neurons;
    double soma_diameter;
    double axon_length, segment_length, sigma_axon_angle; 

    double L; 

    double *X, *Y; // Neuron positions
    double *dendrites_diameters;
    double *axon_lengths;

    int AdjMatrix[NMAX][NMAX]; // Adjacency matrix
    
    // ? #########################################

    // ? ############### INITIALIZATION ###############

    ini_ran(123456789);
    // ini_ran(time(NULL));

    N_neurons = 500;
    L = 2; // Box length (Everything in mm)
    printf("Neurons density: %f neurons/mm^2\n", N_neurons/(L*L));
    
    soma_diameter = 0.015; // Fixed soma diameter (15 microns)
    
    d_mean = 0.3; // Mean of dendrite diameter distribution (Gaussian)
    d_sigma = 0.04; // Standard deviation of dendrite diameter distribution (Gaussian)
    
    l_mean = 1.0; // Mean of axon length distribution (Rayleigh)
    l_sigma = sqrt(2.0/PI); // Sigma parameter for Rayleigh distribution with given mean
    segment_length = 0.01; // Length of each segment in the axon (10 microns)
    sigma_axon_angle = 0.1; // Standard deviation of axon angle deviation (radians)
    
    printf ("Dendrite diameter distribution (Gaussian): mean = %f, sigma = %f\n", d_mean, d_sigma);
    printf ("Axon length distribution (Rayleigh): mean = %f, sigma = %f\n", l_mean, l_sigma);

    X = (double *)malloc(N_neurons * sizeof(double));
    Y = (double *)malloc(N_neurons * sizeof(double));
    dendrites_diameters = (double *)malloc(N_neurons * sizeof(double));
    axon_lengths = (double *)malloc(N_neurons * sizeof(double));

    for(int i = 0; i < N_neurons; i++){
        X[i] = Y[i] = 0.0;
        dendrites_diameters[i] = 0.0;
        axon_lengths[i] = 0.0;
    }

    for(int i = 0; i < NMAX; i++){
        for(int j = 0; j < NMAX; j++){
            AdjMatrix[i][j] = 0;
        }
    }

    // ? #############################################

    // ? ############### FILES ###############

    /* FILE *test_rayleigh;
    if ((test_rayleigh = fopen("test_rayleigh.txt", "w")) == NULL) {
        printf("Error opening file test_rayleigh.txt\n");
        exit(1);
    }

    FILE *test_gaussian;
    if ((test_gaussian = fopen("test_gaussian.txt", "w")) == NULL) {
        printf("Error opening file test_gaussian.txt\n");
        exit(1);
    } */

    FILE *neurons_params;
    if ((neurons_params = fopen("neurons_params.txt", "w")) == NULL) {
        printf("Error opening file neurons_params.txt\n");
        exit(1);
    }
    fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

    /* // ? ############### DISTRIBUTIONS TESTING ###############

    for(int i=0; i<N; i++){
        double u = randomInPR(0.0, 1.0);
        x = inverse_cumulative_rayleigh(u, l_sigma);
        fprintf(test_rayleigh, "%f\n", x);
    }

    for(int i=0; i<N; i++){
        x = d_mean + d_sigma * box_muller();
        fprintf(test_gaussian, "%f\n", x);
    }

    */ // ? ###########################################################


    // ? ############### NETWORK CREATION ###############

    for(int i = 0; i < N_neurons; i++){

        X[i] = randomInPR(-L/2.0, L/2.0);
        Y[i] = randomInPR(-L/2.0, L/2.0);

        dendrites_diameters[i] = d_mean + d_sigma * box_muller();
        double u = randomInPR(0.0, 1.0);
        axon_lengths[i] = inverse_cumulative_rayleigh(u, l_sigma);
        fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
    
    }

    for(int i = 0; i < N_neurons; i++){}


    fclose(neurons_params);
    return 0;
}