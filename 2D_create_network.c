#include "2D_head.h"

int main(int argc, char *argv[]) {

    if (argc < 2){
        printf("Usage: %s alpha\n", argv[0]);
        return 1;
    }

    
    ini_ran(time(NULL));
    
    double *X, *Y;
    
    double L, rho; int N_neurons;
    double soma_diameter;
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh;
    double segment_length, sigma_axon_angle;
    int n_centers;
    double base_sigma;
    double alpha;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];
    char filename_simulation[MAX_STRING_LENGTH];

    FILE *axon_simulation;
    FILE *neurons_params;
    
// !! FILENAME TO LOAD PARAMETERS AND NEURON CONFIGURATIONS

    sprintf(filename, "2D_initial_configurations/2D_neurons_params_random_L2.0_rho125_l1.00.txt");

// !! FILENAME TO LOAD PARAMETERS AND NEURON CONFIGURATIONS

    if ((neurons_params = fopen(filename, "r")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    fscanf(neurons_params, "L = %lf\n", &L);
    fscanf(neurons_params, "rho = %lf\n", &rho);
    fscanf(neurons_params, "soma_diameter = %lf\n", &soma_diameter);
    fscanf(neurons_params, "d_mean = %lf\n", &d_mean);
    fscanf(neurons_params, "d_sigma = %lf\n", &d_sigma);
    fscanf(neurons_params, "l_mean = %lf\n", &l_mean);
    fscanf(neurons_params, "l_sigma = %lf\n", &l_sigma);
    fscanf(neurons_params, "sigma_rayleigh = %lf\n", &sigma_rayleigh);
    fscanf(neurons_params, "segment_length = %lf\n", &segment_length);
    fscanf(neurons_params, "sigma_axon_angle = %lf\n", &sigma_axon_angle);
    fscanf(neurons_params, "nc = %d\n", &n_centers);
    fscanf(neurons_params, "base_sigma = %lf\n", &base_sigma);

    N_neurons = (int)(rho * L * L);

    printf("Parameters loaded from file: %s\n", filename);
    printf("L = %.1lf\n", L);
    printf("rho = %.0f --> N = %d\n", rho, N_neurons);
    printf("soma_diameter = %.3lf\n", soma_diameter);
    printf("d_mean = %.3lf\n", d_mean);
    printf("d_sigma = %.3lf\n", d_sigma);
    printf("l_mean = %.3lf\n", l_mean);
    printf("l_sigma = %.3lf\n", l_sigma);
    printf("sigma_rayleigh = %.3lf\n", sigma_rayleigh);
    printf("segment_length = %.3lf\n", segment_length);
    printf("sigma_axon_angle = %.3lf\n", sigma_axon_angle);
    printf("nc = %d\n", n_centers);
    printf("base_sigma = %.3lf\n", base_sigma);

    fgets(filename, MAX_STRING_LENGTH, neurons_params); // Read the header line
    printf("Heades line skipped: %s\n", filename);

    X = (double *)malloc(N_neurons * sizeof(double));
    Y = (double *)malloc(N_neurons * sizeof(double));

    somas = (double *)malloc(N_neurons * sizeof(double));
    dendrites_diameters = (double *)malloc(N_neurons * sizeof(double));
    axon_lengths = (double *)malloc(N_neurons * sizeof(double));

    for(int i = 0; i < N_neurons; i++){
        fscanf(neurons_params, "%lf\t%lf\t%lf\t%lf\t%lf\n", &X[i], &Y[i], &somas[i], &dendrites_diameters[i], &axon_lengths[i]);
    }

    fclose(neurons_params);

    // alpha = 0.2; 
    alpha = atof(argv[1]); // ! CONNECTION PROBABILITY

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

    if (n_centers == 0){
        sprintf(filename_simulation, "2D_axon_simulation/2D_axon_positions_random_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", L, rho, l_mean, alpha);
    }
    else{
        sprintf(filename_simulation, "2D_axon_simulation/2D_axon_positions_agg_nc%d_s%.2f_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", n_centers, base_sigma, L, rho, l_mean, alpha);
    }

    if ((axon_simulation = fopen(filename_simulation, "w")) == NULL) {
        printf("Error opening file %s\n", filename_simulation);
        exit(1);
    }

    fprintf(axon_simulation, "Neuron_Index\tSegment_Index\tStart_X\tStart_Y\tEnd_X\tEnd_Y\n");

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK

    if (n_centers == 0){
        sprintf(filename, "2D_created_networks/2D_adjacency_matrix_random_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", L, rho, l_mean, alpha);
    }
    else{
        sprintf(filename, "2D_created_networks/2D_adjacency_matrix_agg_nc%d_s%.2f_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", n_centers, base_sigma, L, rho, l_mean, alpha);
    }

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK 

    // ? CREATION OF THE NETWORK

    int trials, links;

    trials = 0; 
    links = 0;

    uint8_t *AdjMatrix_flat = (uint8_t *)malloc(N_neurons * N_neurons * sizeof(uint8_t));
    memset(AdjMatrix_flat, 0, N_neurons * N_neurons * sizeof(uint8_t)); // Initialize adjacency matrix with zeros

    for(int i = 0; i < N_neurons; i++){
        printf("Creating connections for neuron %d/%d\r", i+1, N_neurons);
        double initial_segment_vector_x, initial_segment_vector_y;
        double initial_angle = randomInPR(0.0, 2*PI); // Initial angle of the axon (randomly distributed between 0 and 2*pi)

        initial_segment_vector_x = X[i] + (soma_diameter/2.0) * cos(initial_angle);
        initial_segment_vector_y = Y[i] + (soma_diameter/2.0) * sin(initial_angle);
        PBC(&initial_segment_vector_x, L);
        PBC(&initial_segment_vector_y, L);

        int number_of_segments = (int)(axon_lengths[i] / segment_length);

        for(int k = 0; k < number_of_segments; k++){

            double dx, dy; 
            double angle; 
            double end_segment_vector_x, end_segment_vector_y;

            if (k ==0){
                angle = initial_angle;
            }
            else{
                angle = initial_angle + sigma_axon_angle * box_muller(); // New angle of the axon segment (Gaussian deviation from the initial angle)
            }

            new_vector_segment(segment_length, angle, &dx, &dy); // Calculate the vector of the new segment based on the angle and segment length
            initial_angle = angle; // Update the initial angle for the next segment

            for(int j = 0; j<N_neurons; j++){
                if (i != j){
                    if(AdjMatrix_flat[i*N_neurons + j] == 0){ // If there is no connection yet
                        if (new_axon_intersection(X, Y, dendrites_diameters, i, j, &initial_segment_vector_x, &initial_segment_vector_y, L, dx, dy)){ // Check if the new segment intersects with the dendrites of neuron j
                            trials ++;
                            if (randomInPR(0.0, 1.0) < alpha){ // Connection probability alpha
                                AdjMatrix_flat[i*N_neurons + j] = 1; // Create connection
                                links ++;
                            }
                        }
                    }
                }
            }

            end_segment_vector_x = initial_segment_vector_x + dx;
            end_segment_vector_y = initial_segment_vector_y + dy;

            PBC(&end_segment_vector_x, L);
            PBC(&end_segment_vector_y, L);

            fprintf(axon_simulation, "%d\t%d\t%f\t%f\t%f\t%f\n", i, k, initial_segment_vector_x, initial_segment_vector_y, end_segment_vector_x, end_segment_vector_y);

            initial_segment_vector_x = end_segment_vector_x;
            initial_segment_vector_y = end_segment_vector_y;
        }
    }

    fclose(axon_simulation);

    printf("\nTrials (intersections) = %d\n", trials);
    printf("Links created          = %d\n", links);
    printf("Acceptance ratio       = %.6f\n", (double)links / (double)trials);

    FILE *adj_matrix_file;
    if ((adj_matrix_file = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    for(int i = 0; i < N_neurons; i++){
        for(int j = 0; j < N_neurons; j++){
            fprintf(adj_matrix_file, "%d ", AdjMatrix_flat[i*N_neurons + j]);
        }
        fprintf(adj_matrix_file, "\n");
    }
    
    fclose(adj_matrix_file);

    free(X);
    free(Y);
    free(somas);
    free(dendrites_diameters);
    free(axon_lengths);
    free(AdjMatrix_flat);

    return 0;
}

// double P = 0.0;
// for(int j = 0; j<N_neurons; j++){
//     if (i != j){
//         if(AdjMatrix_flat[i*N_neurons + j] == 0){
//             double dxj = min_image(end_segment_x - X[j], L);
//             double dyj = min_image(end_segment_y - Y[j], L);
//             P = dxj*dxj + dyj*dyj;
//             if (P <= (dendrites_diameters[j]/2.0)*(dendrites_diameters[j]/2.0)){
//                 double omega = randomInPR(0.0, 1.0);
//                 if (omega < alpha){ // Connection probability alpha
//                     AdjMatrix_flat[i*N_neurons + j] = 1;
//                 }
//             }
//         }
//     }
// }