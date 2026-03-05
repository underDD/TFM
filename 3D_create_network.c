#include "3D_head.h"

int main(int argc, char *argv[]) {

    if (argc < 2){
        printf("Usage: %s alpha\n", argv[0]);
        return 1;
    }

    ini_ran(time(NULL));

    double *X, *Y, *Z;

    double L, rho; int N;
    double soma_diameter;
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh;
    double segment_length, sigma_pol, sigma_azi;
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

    sprintf(filename, "3D_initial_configurations/3D_neurons_params_random_L2.0_rho125_l1.00.txt");

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
    fscanf(neurons_params, "sigma_pol = %lf\n", &sigma_pol);
    fscanf(neurons_params, "sigma_azi = %lf\n", &sigma_azi);
    fscanf(neurons_params, "nc = %d\n", &n_centers);
    fscanf(neurons_params, "base_sigma = %lf\n", &base_sigma);
    
    N = (int)(rho * L * L * L);

    printf("Parameters loaded from file: %s\n", filename);
    printf("L = %.1lf\n", L);
    printf("rho = %.0f --> N = %d\n", rho, N);
    printf("soma_diameter = %.3lf\n", soma_diameter);
    printf("d_mean = %.3lf\n", d_mean);
    printf("d_sigma = %.3lf\n", d_sigma);
    printf("l_mean = %.3lf\n", l_mean);
    printf("l_sigma = %.3lf\n", l_sigma);
    printf("sigma_rayleigh = %.3lf\n", sigma_rayleigh);
    printf("segment_length = %.3lf\n", segment_length);
    printf("sigma_pol = %.3lf\n", sigma_pol);
    printf("sigma_azi = %.3lf\n", sigma_azi);
    printf("nc = %d\n", n_centers);
    printf("base_sigma = %.3lf\n", base_sigma);

    fgets(filename, MAX_STRING_LENGTH, neurons_params); // Skip header line
    printf("Header line skipped: %s\n", filename);

    X = (double *)malloc(N * sizeof(double));
    Y = (double *)malloc(N * sizeof(double));
    Z = (double *)malloc(N * sizeof(double));

    somas = (double *)malloc(N * sizeof(double));
    dendrites_diameters = (double *)malloc(N * sizeof(double));
    axon_lengths = (double *)malloc(N * sizeof(double));

    for(int i = 0; i < N; i++){
        fscanf(neurons_params, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &X[i], &Y[i], &Z[i], &somas[i], &dendrites_diameters[i], &axon_lengths[i]);
    }

    fclose(neurons_params);

    // alpha = 0.2; 
    alpha = atof(argv[1]); // ! CONNECTION PROBABILITY

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

    if (n_centers == 0){
        sprintf(filename_simulation, "3D_axon_simulation/3D_axon_positions_random_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", L, rho, l_mean, alpha);
    }
    else{
        sprintf(filename_simulation, "3D_axon_simulation/3D_axon_positions_agg_nc%d_s%.2f_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", n_centers, base_sigma, L, rho, l_mean, alpha);
    }

    if ((axon_simulation = fopen(filename_simulation, "w")) == NULL) {
        printf("Error opening file %s\n", filename_simulation);
        exit(1);
    }

    fprintf(axon_simulation, "Neuron_Index\tSegment_Index\tStart_X\tStart_Y\tStart_Z\tEnd_X\tEnd_Y\tEnd_Z\n");

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK

    if (n_centers == 0){
        sprintf(filename, "3D_created_networks/3D_adjacency_matrix_random_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", L, rho, l_mean, alpha);
    }
    else{
        sprintf(filename, "3D_created_networks/3D_adjacency_matrix_agg_nc%d_s%.2f_L%.1lf_rho%.0f_l%.2f_a%.3f.txt", n_centers, base_sigma, L, rho, l_mean, alpha);
    }

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK

    // ? CREATION OF THE NETWORK
    
    uint8_t *AdjMatrix_flat = (uint8_t *)malloc(N * N * sizeof(uint8_t));
    memset(AdjMatrix_flat, 0, N * N * sizeof(uint8_t));

    int trials, links;

    trials = 0; 
    links = 0;

    for(int i = 0; i < N; i++){
        printf("Creating connections for neuron %d/%d\r", i+1, N);
        double initial_segment_vector_x, initial_segment_vector_y, initial_segment_vector_z;
        double initial_angle_azi = randomInPR(0.0, 2.0*PI);
        double initial_angle_pol = acos(randomInPR(-1.0, 1.0)); // Polar angle between 0 and pi

        initial_segment_vector_x = X[i] + (soma_diameter/2.0) * cos(initial_angle_azi) * sin(initial_angle_pol);
        initial_segment_vector_y = Y[i] + (soma_diameter/2.0) * sin(initial_angle_azi) * sin(initial_angle_pol);
        initial_segment_vector_z = Z[i] + (soma_diameter/2.0) * cos(initial_angle_pol);
        PBC(&initial_segment_vector_x, L);
        PBC(&initial_segment_vector_y, L);
        PBC(&initial_segment_vector_z, L);

        int number_of_segments = (int)(axon_lengths[i] / segment_length);

        for(int k = 0; k < number_of_segments; k++){

            double dx, dy, dz; 
            double angle_pol, angle_azi;
            double end_segment_vector_x, end_segment_vector_y, end_segment_vector_z;

            if (k == 0){
                angle_azi = initial_angle_azi;
                angle_pol = initial_angle_pol;
            }
            else{
                angle_azi = initial_angle_azi + box_muller() * sigma_azi; // Standard deviation of sigma_azi radians
                angle_pol = initial_angle_pol + box_muller() * sigma_pol; // Standard deviation of sigma_pol radians

                angle_pol = fmod(angle_pol, 2.0*PI);
                if (angle_pol < 0.0) angle_pol += 2.0*PI;
                if (angle_pol > PI) angle_pol = 2.0*PI - angle_pol; // Keep polar angle between 0 and pi

                angle_azi = fmod(angle_azi, 2.0*PI);
                if (angle_azi < 0.0) angle_azi += 2.0*PI; // Keep azimuthal angle between 0 and 2pi

            }

            new_vector_segment_3D(segment_length, angle_azi, angle_pol, &dx, &dy, &dz);
            initial_angle_azi = angle_azi;
            initial_angle_pol = angle_pol;

            for(int j = 0; j < N; j++){
                if (i != j){
                    if(AdjMatrix_flat[i*N + j] == 0){
                        if (new_axon_intersection_3D(X, Y, Z, dendrites_diameters, i, j, &initial_segment_vector_x, &initial_segment_vector_y, &initial_segment_vector_z, L, dx, dy, dz)){
                            trials ++;
                            if (randomInPR(0.0, 1.0) < alpha){ // Connection probability alpha
                                AdjMatrix_flat[i*N + j] = 1;
                                links ++;
                            }
                        }
                    }
                }
            }

            end_segment_vector_x = initial_segment_vector_x + dx;
            end_segment_vector_y = initial_segment_vector_y + dy;
            end_segment_vector_z = initial_segment_vector_z + dz;

            PBC(&end_segment_vector_x, L);
            PBC(&end_segment_vector_y, L);
            PBC(&end_segment_vector_z, L);

            fprintf(axon_simulation, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, k, initial_segment_vector_x, initial_segment_vector_y, initial_segment_vector_z, end_segment_vector_x, end_segment_vector_y, end_segment_vector_z);

            initial_segment_vector_x = end_segment_vector_x;
            initial_segment_vector_y = end_segment_vector_y;
            initial_segment_vector_z = end_segment_vector_z;

        }
    }

    fclose(axon_simulation);

    printf("\nTrials (intersections) = %d\n", trials);
    printf("Links created          = %d\n", links);
    printf("Acceptance ratio       = %.6f\n", (double)links / (double)trials);

    // Save adjacency matrix to file
    FILE *adj_matrix_file;
    if ((adj_matrix_file = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            fprintf(adj_matrix_file, "%d ", AdjMatrix_flat[i*N + j]);
        }
        fprintf(adj_matrix_file, "\n");
    }
    fclose(adj_matrix_file);


    free(X);
    free(Y);
    free(Z);
    free(somas);
    free(dendrites_diameters);
    free(axon_lengths);
    free(AdjMatrix_flat);
    
    return 0;
}

/*
    double P = 0.0;
    for(int j = 0; j<N_neurons; j++){
        if (i != j){
            if(AdjMatrix_flat[i*N_neurons + j] == 0){
                double dxj = min_image(end_segment_x - X[j], L);
                double dyj = min_image(end_segment_y - Y[j], L);
                double dzj = min_image(end_segment_z - Z[j], L);
                P = dxj*dxj + dyj*dyj + dzj*dzj;
                if (P <= (dendrites_diameters[j]/2.0)*(dendrites_diameters[j]/2.0)){
                    double omega = randomInPR(0.0, 1.0);
                    if (omega < alpha){ // Connection probability alpha
                        AdjMatrix_flat[i*N_neurons + j] = 1;
                    }
                }
            }
        }
    }
*/
