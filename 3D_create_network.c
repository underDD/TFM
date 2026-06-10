#include "3D_head.h"

int main(int argc, char *argv[]) {

    if (argc < 16) {
        printf("Usage: %s filename L rho soma_diameter d_mean d_sigma l_mean segment_length sigma_pol sigma_azi n_centers base_sigma BC_type geometry seed\n", argv[0]);
        return 1;
    }

    // ini_ran(seed);

    double *X, *Y, *Z;

    double L, rho; long N;
    double soma_diameter;
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh;
    double segment_length, sigma_pol, sigma_azi;
    int n_centers;
    double base_sigma;
    double alpha; 
    double R, height;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];
    char filename_simulation[MAX_STRING_LENGTH];
    char BC_type[MAX_STRING_LENGTH];
    char geometry[MAX_STRING_LENGTH];

    // strcpy(BC_type, "PBC"); // Boundary conditions type (PBC: Periodic Boundary Conditions, RBC: Reflective Boundary Conditions)

    FILE *axon_simulation;
    FILE *neurons_params;
    
// !! FILENAME TO LOAD NEURON CONFIGURATION

    // sprintf(filename, "3D_initial_configurations/3D_neurons_params_random_L2.0_rho125_l1.00.txt");
    strcpy(filename, argv[1]);
    L = atof(argv[2]);
    rho = atof(argv[3]);
    soma_diameter = atof(argv[4]);
    d_mean = atof(argv[5]);
    d_sigma = atof(argv[6]);
    l_mean = atof(argv[7]);
    segment_length = atof(argv[8]);
    sigma_pol = atof(argv[9]);
    sigma_azi = atof(argv[10]);
    n_centers = atoi(argv[11]);
    base_sigma = atof(argv[12]);
    strcpy(BC_type, argv[13]);
    strcpy(geometry, argv[14]);
    int seed = atoi(argv[15]);

    if (seed!= 0){
        ini_ran(seed);
    }
    else{
        ini_ran((int)time(NULL));
    }

// !! FILENAME TO LOAD PARAMETERS AND NEURON CONFIGURATIONS

    if ((neurons_params = fopen(filename, "r")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    
    if (strcmp(geometry, "cube") == 0 || strcmp(geometry, "cubeBiased") == 0 || strcmp(geometry, "cubeLongRange") == 0){
        N = (long)(rho * L * L * L);
    }
    else if (strcmp(geometry, "cylinder") == 0){
        R = L / 2; // Calculate the radius of the circle to maintain the same area as the square
        height = L; // Keep the same height as the side of the square
        N = (long)(rho * PI * R * R * height);
    }

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
    printf("base_sigma = %.4lf\n", base_sigma);
    printf("BC_type = %s\n", BC_type);
    printf("Geometry = %s\n", geometry);

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
    // alpha = atof(argv[1]); // ! CONNECTION PROBABILITY

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

    if (n_centers == 0){
        sprintf(filename_simulation, "3D_axon_simulation/3D_axon_positions_%s_%s_random_L%.1lf_rho%.0f_l%.2f_d%.2f.txt", geometry, BC_type, L, rho, l_mean, d_mean);
    }
    else{
        sprintf(filename_simulation, "3D_axon_simulation/3D_axon_positions_%s_%s_agg_nc%d_s%.4f_L%.1lf_rho%.0f_l%.2f_d%.2f.txt", geometry, BC_type, n_centers, base_sigma, L, rho, l_mean, d_mean);
    }

    if ((axon_simulation = fopen(filename_simulation, "w")) == NULL) {
        printf("Error opening file %s\n", filename_simulation);
        exit(1);
    }

    if (strcmp(BC_type, "PBC") == 0){
        fprintf(axon_simulation, "Neuron_Index\tSegment_Index\tStart_X\tStart_Y\tStart_Z\tEnd_X\tEnd_Y\tEnd_Z\n");
    }
    else{
        fprintf(axon_simulation, "Neuron_Index\tSegment_Index\tSubsegment_Index\tStart_X\tStart_Y\tStart_Z\tEnd_X\tEnd_Y\tEnd_Z\n");
    }

// ! FILENAME TO SAVE THE POSITIONS OF THE AXONS TO THE SIMULATION

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK

    if (n_centers == 0){
        sprintf(filename, "3D_created_networks/3D_adjacency_matrix_%s_%s_random_L%.1lf_rho%.0f_l%.2f_d%.2f.txt", geometry, BC_type, L, rho, l_mean, d_mean);
    }
    else{
        sprintf(filename, "3D_created_networks/3D_adjacency_matrix_%s_%s_agg_nc%d_s%.4f_L%.1lf_rho%.0f_l%.2f_d%.2f.txt", geometry, BC_type, n_centers, base_sigma, L, rho, l_mean, d_mean);
    }

// ! FILENAME TO SAVE THE ADJACENCY MATRIX OF THE CREATED NETWORK

    // ? CREATION OF THE NETWORK
    
    uint8_t *AdjMatrix_flat = (uint8_t *)malloc(N * N * sizeof(uint8_t));
    memset(AdjMatrix_flat, 0, N * N * sizeof(uint8_t));

    int trials, links;

    trials = 0; 
    links = 0;

    for(int i = 0; i < N; i++){
        // printf("Creating connections for neuron %d/%d\r", i+1, N);
        
        if (seed!= 0){
            ini_ran(seed + i*1000); // Re-seed the random number generator for each neuron to ensure different growth patterns
        }
        
        double xmin = -L/2.0, xmax = L/2.0;
        double ymin = -L/2.0, ymax = L/2.0;
        double zmin = -L/2.0, zmax = L/2.0;

        fflush(stdout);
        double initial_segment_vector_x, initial_segment_vector_y, initial_segment_vector_z;
        double initial_angle_azi = randomInPR(0.0, 2.0*PI);
        double initial_angle_pol;
        double sigma_pol_init_biased = 0.01;
        double p_vertical_axons = 0.5;  // 10% pueden crecer más en Z, el resto crece más en el plano XY
        int vertical_axon = 0;

        if (strcmp(geometry, "cubeBiased") == 0){
            if (randomInPR(0.0, 1.0) < p_vertical_axons){
                vertical_axon = 1;
            }
        }

        if (strcmp(geometry, "cubeBiased") == 0){
            if (vertical_axon){
                initial_angle_pol = acos(randomInPR(-1.0, 1.0));
            }
            else{
                initial_angle_pol = PI / 2.0 + box_muller() * sigma_pol_init_biased;
            }
        }
        else{
            initial_angle_pol = acos(randomInPR(-1.0, 1.0));
        }

        double left_band_width = 0.30;   // grosor del lateral izquierdo
        double p_long_range = 0.30;      // fracción de neuronas del lateral que lanzan axón largo
        int long_range_axon = 0;

        if (strcmp(geometry, "cubeLongRange") == 0) {

            if (X[i] < xmin + left_band_width) {
                if (randomInPR(0.0, 1.0) < p_long_range) {
                    long_range_axon = 1;
                }
            }
        }

        if (long_range_axon) {
            initial_angle_azi = 0.0;       // dirección +X
            initial_angle_pol = PI / 2.0;  // plano XY
        }
        else {
            initial_angle_azi = randomInPR(0.0, 2.0*PI);

            if (strcmp(geometry, "cubeBiased") == 0) {

                if (vertical_axon) {
                    initial_angle_pol = acos(randomInPR(-1.0, 1.0));
                }
                else {
                    initial_angle_pol = PI / 2.0 + box_muller() * sigma_pol_init_biased;
                }
            }
            else {
                initial_angle_pol = acos(randomInPR(-1.0, 1.0));
            }
        }



        initial_segment_vector_x = X[i] + (soma_diameter/2.0) * cos(initial_angle_azi) * sin(initial_angle_pol);
        initial_segment_vector_y = Y[i] + (soma_diameter/2.0) * sin(initial_angle_azi) * sin(initial_angle_pol);
        initial_segment_vector_z = Z[i] + (soma_diameter/2.0) * cos(initial_angle_pol);
        
        if (strcmp(BC_type, "PBC") == 0){
            PBC(&initial_segment_vector_x, L);
            PBC(&initial_segment_vector_y, L);
            PBC(&initial_segment_vector_z, L);
        }
        else{
            if (initial_segment_vector_x < xmin) initial_segment_vector_x = xmin;
            if (initial_segment_vector_x > xmax) initial_segment_vector_x = xmax;
            if (initial_segment_vector_y < ymin) initial_segment_vector_y = ymin;
            if (initial_segment_vector_y > ymax) initial_segment_vector_y = ymax;
            if (initial_segment_vector_z < zmin) initial_segment_vector_z = zmin;
            if (initial_segment_vector_z > zmax) initial_segment_vector_z = zmax;
        }

        int number_of_segments = (int)(axon_lengths[i] / segment_length);

        for(int k = 0; k < number_of_segments; k++){

            double dx, dy, dz; 
            double angle_pol, angle_azi;
            double end_segment_vector_x, end_segment_vector_y, end_segment_vector_z;

            if (long_range_axon) {
                angle_azi = 0.0;
                angle_pol = PI / 2.0;
            }
            else if (k == 0) {
                angle_azi = initial_angle_azi;
                angle_pol = initial_angle_pol;
            }
            else{
                angle_azi = initial_angle_azi + box_muller() * sigma_azi; // Standard deviation of sigma_azi radians
                if (strcmp(geometry, "cubeBiased") == 0){

                    if (vertical_axon){
                        angle_pol = initial_angle_pol + box_muller() * sigma_pol;
                    }
                    else{
                        double sigma_pol_biased = 0.01;
                        double lambda_plane = 0.5;

                        angle_pol = initial_angle_pol
                                + box_muller() * sigma_pol_biased
                                - lambda_plane * (initial_angle_pol - PI / 2.0);
                    }
                }
                else{
                    angle_pol = initial_angle_pol + box_muller() * sigma_pol;
                }
                angle_pol = fmod(angle_pol, 2.0*PI);
                if (angle_pol < 0.0) angle_pol += 2.0*PI;
                if (angle_pol > PI) angle_pol = 2.0*PI - angle_pol; // Keep polar angle between 0 and pi

                angle_azi = fmod(angle_azi, 2.0*PI);
                if (angle_azi < 0.0) angle_azi += 2.0*PI; // Keep azimuthal angle between 0 and 2pi

            }

            new_vector_segment_3D(segment_length, angle_azi, angle_pol, &dx, &dy, &dz);
            initial_angle_azi = angle_azi;
            initial_angle_pol = angle_pol;

            if (strcmp(BC_type, "PBC") == 0){
                for(int j = 0; j < N; j++){
                    if (i != j){
                    
                            if (new_axon_intersection_3D(X, Y, Z, dendrites_diameters, i, j, &initial_segment_vector_x, &initial_segment_vector_y, &initial_segment_vector_z, L, dx, dy, dz)){
                                trials ++;
                                 // Connection probability alpha
                                AdjMatrix_flat[i*N + j] += 1;
                                links ++;
                                
                            }
                    }
                }
                end_segment_vector_x = initial_segment_vector_x + dx;
                end_segment_vector_y = initial_segment_vector_y + dy;
                end_segment_vector_z = initial_segment_vector_z + dz;
    
                PBC(&end_segment_vector_x, L);
                PBC(&end_segment_vector_y, L);
                PBC(&end_segment_vector_z, L);
            }
            else if ((strcmp(geometry, "cube") == 0 ||
                    strcmp(geometry, "cubeBiased") == 0 ||
                    strcmp(geometry, "cubeLongRange") == 0)
                    && strcmp(BC_type, "SW") == 0){
                sticky_walls3D(initial_segment_vector_x, initial_segment_vector_y, initial_segment_vector_z, 
                                &dx, &dy, &dz, L, &end_segment_vector_x, &end_segment_vector_y, &end_segment_vector_z, X, Y, Z, 
                                dendrites_diameters, i, N, AdjMatrix_flat, &trials, &links, axon_simulation, i, k);
            
                if (dx*dx + dy*dy + dz*dz > 1e-15){
                    initial_angle_azi = atan2(dy, dx);
                    initial_angle_pol = acos(dz / sqrt(dx*dx + dy*dy + dz*dz));
                }
            }
            else if (strcmp(geometry, "cylinder") == 0 && strcmp(BC_type, "SW") == 0){
                int stop_axon = sticky_cylinder3D(
                    initial_segment_vector_x,
                    initial_segment_vector_y,
                    initial_segment_vector_z,
                    &dx, &dy, &dz,
                    L,
                    &end_segment_vector_x,
                    &end_segment_vector_y,
                    &end_segment_vector_z,
                    X, Y, Z,
                    dendrites_diameters,
                    i, N,
                    AdjMatrix_flat,
                    &trials, &links,
                    axon_simulation,
                    i, k
                );

                if (dx*dx + dy*dy + dz*dz > 1e-15) {
                    initial_angle_azi = atan2(dy, dx);
                    initial_angle_pol = acos(dz / sqrt(dx*dx + dy*dy + dz*dz));
                }

                if (stop_axon) {
                    initial_segment_vector_x = end_segment_vector_x;
                    initial_segment_vector_y = end_segment_vector_y;
                    initial_segment_vector_z = end_segment_vector_z;
                    break;
                }
            
            }

            // fprintf(axon_simulation, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, k, initial_segment_vector_x, initial_segment_vector_y, initial_segment_vector_z, end_segment_vector_x, end_segment_vector_y, end_segment_vector_z);

            initial_segment_vector_x = end_segment_vector_x;
            initial_segment_vector_y = end_segment_vector_y;
            initial_segment_vector_z = end_segment_vector_z;

        }
    }
    printf("\n");
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
