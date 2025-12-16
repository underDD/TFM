#include "head.h"

int main(){
    printf("patata\n");

    // ? ############### VARIABLES ###############

    double l_mean, l_sigma; 
    double d_mean, d_sigma;

    int N_neurons;
    double soma_diameter;
    double segment_length, sigma_axon_angle; 
    double alpha;

    double L; 

    double *X, *Y; // Neuron positions
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];

    int **AdjMatrix = malloc(NMAX * sizeof(int *));
    if (AdjMatrix == NULL) {
        perror("Error allocating memory for AdjMatrix");
        exit(1); // Maneja el error si la asignación falla
    }

    for (int i = 0; i < NMAX; i++) {
        AdjMatrix[i] = malloc(NMAX * sizeof(int));
        if (AdjMatrix[i] == NULL) {
            perror("Error allocating memory for AdjMatrix rows");
            exit(1);
        }
    }

    
    // ? #########################################

    // ? ############### INITIALIZATION ###############

    ini_ran(123456789);
    // ini_ran(time(NULL));

    N_neurons = 500;
    L = 2; // Box length (Everything in mm)
    printf("Neurons density: %f neurons/mm^2\n", N_neurons/(L*L));
    
    soma_diameter = 0.015; // Fixed soma diameter (15 microns)
    
    d_mean = 0.3; // Mean of dendrite diameter distribution (Gaussian)
    d_sigma = sqrt(0.04); // Standard deviation of dendrite diameter distribution (Gaussian)
    
    double sigma_rayleigh;
    l_mean = 1.0;
    l_sigma = sqrt(0.8); // Sigma parameter for Rayleigh distribution with given mean
    // sigma_rayleigh = l_mean*sqrt(2/PI); // Mean of axon length distribution (Rayleigh)
    sigma_rayleigh = l_sigma*sqrt(2/(4 - PI)); // Recalculate sigma from desired standard deviation
    l_mean = sigma_rayleigh*sqrt(PI/2.0); // Calculate mean from sigma
    segment_length = 0.01; // Length of each segment in the axon (10 microns)
    sigma_axon_angle = sqrt(0.1); // Standard deviation of axon angle deviation (radians)

    alpha = 0.01;
    
    printf("Dendrite diameter distribution (Gaussian): mean = %f, sigma = %f\n", d_mean, d_sigma);
    printf("Axon length distribution (Rayleigh): mean = %f, sigma = %f\n", l_mean, l_sigma);

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

    // FILE *test_rayleigh;
    // if ((test_rayleigh = fopen("test_rayleigh.txt", "w")) == NULL) {
    //     printf("Error opening file test_rayleigh.txt\n");
    //     exit(1);
    // }

    // FILE *test_gaussian;
    // if ((test_gaussian = fopen("test_gaussian.txt", "w")) == NULL) {
    //     printf("Error opening file test_gaussian.txt\n");
    //     exit(1);
    // } 

    sprintf(filename, "networks/neurons_params_a%.2f.txt", alpha);
    FILE *neurons_params;
    if ((neurons_params = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

    sprintf(filename, "networks/adjacency_matrix_a%.2f.txt", alpha);

    FILE *adjacency_matrix;
    if ((adjacency_matrix = fopen(filename, "w")) == NULL){
        printf("Error opening file %s\n", filename);
        exit(1);
    }

     // ? ############### DISTRIBUTIONS TESTING ###############

    // int N = 10000;
    // double x;

    // for(int i=0; i<N; i++){
    //     double u = randomInPR(0.0, 1.0);
    //     x = inverse_cumulative_rayleigh(u, sigma_rayleigh);
    //     fprintf(test_rayleigh, "%f\n", x);
    // }

    // for(int i=0; i<N; i++){
    //     x = d_mean + d_sigma * box_muller();
    //     fprintf(test_gaussian, "%f\n", x);
    // }

     // ? ###########################################################


    // ? ############### NETWORK CREATION ###############

    bool overlap;
    for(int i = 0; i < N_neurons; i++){ // Only to place the neurons in the plane, without making links yet
        
        do{
            overlap = false;
            X[i] = randomInPR(-L/2.0, L/2.0);
            Y[i] = randomInPR(-L/2.0, L/2.0);
            for(int j = 0; j < i; j++){
                double dist = sqrt_distance(X[i]-X[j], Y[i]-Y[j]);
                if (dist < soma_diameter){
                    overlap = true;
                    break;
                }
            }
        }while(overlap);

        dendrites_diameters[i] = d_mean + d_sigma * box_muller();
        double u = randomInPR(0.0, 1.0);
        axon_lengths[i] = inverse_cumulative_rayleigh(u, l_sigma);
        fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
    
    }

    printf("Neurons placed. Creating connections...\n");

    for(int i = 0; i < N_neurons; i++){
        double segment_vector_x, segment_vector_y;
        double initial_angle = randomInPR(-PI, PI);
        segment_vector_x = X[i] + (soma_diameter/2.0) * cos(initial_angle);
        segment_vector_y = Y[i] + (soma_diameter/2.0) * sin(initial_angle);
        int number_of_segments = (int)(axon_lengths[i] / segment_length);
        for(int j = 0; j<N_neurons; j++){
            if (i != j){
                for(int k = 0; k < number_of_segments; k++){
                    bool intersection = new_axon_intersection(X, Y, dendrites_diameters, i, j, &segment_vector_x, &segment_vector_y, segment_length, sigma_axon_angle, L);
                    if (intersection){
                        if (randomInPR(0.0, 1.0) < alpha){ // Connection probability alpha
                            AdjMatrix[i][j] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < N_neurons; i++){
        for(int j = 0; j < N_neurons; j++){
            fprintf(adjacency_matrix, "%d\t", AdjMatrix[i][j]);
        }
        fprintf(adjacency_matrix, "\n");
    }

    // ? #############################################

    // fclose(test_rayleigh);
    // fclose(test_gaussian); 

    // Al final, libera la memoria
    for (int i = 0; i < NMAX; i++) {
        free(AdjMatrix[i]);
    }
    free(AdjMatrix);
    fclose(neurons_params);

    free(X);
    free(Y);
    free(dendrites_diameters);
    free(axon_lengths);

    return 0;
}