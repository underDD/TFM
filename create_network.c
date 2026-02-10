#include "head.h"

int main(){

    // ? ############### VARIABLES ###############

    double l_mean, l_sigma; 
    double d_mean, d_sigma;

    int N_neurons;
    double soma_diameter;
    double segment_length, sigma_axon_angle; 

    double L; 

    double *X, *Y; // Neuron positions
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];

    bool aggregation;

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

    // ini_ran(123456789);
    ini_ran(time(NULL));

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

    double alpha_0;
    double alpha_f;
    double d_alpha;
    double alpha;

    alpha_0 = 0.001;
    alpha_f = 1.0;
    d_alpha = 0.001;

    aggregation = true; // true for aggregation, false for no aggregation

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

    int number_of_intersections = 0;
    double mean_intersections = 0.0;

    int runs = 20;

    for(int run = 10; run < runs; run++){

        ini_ran(time(NULL) + run*1000);
        printf("Run %d/%d\n", run+1, runs);

        sprintf(filename, "networks/adjacency_matrix_ng1_run%03d.bin", run);
        FILE *fb = fopen(filename, "wb");
        if (!fb) { perror("Error opening bin file"); exit(1); }

        int32_t N = N_neurons; // int32_t - integer with sign of 32 bites
        int32_t M = (int32_t)((alpha_f - alpha_0)/d_alpha) + 1;
        fwrite(&N, sizeof(int32_t), 1, fb);
        fwrite(&M, sizeof(int32_t), 1, fb);
        fwrite(&alpha_0, sizeof(double), 1, fb);
        fwrite(&d_alpha, sizeof(double), 1, fb);

        uint8_t *AdjMatrix_flat; // uint8_t - integer without sign of 8 bites
        AdjMatrix_flat = (uint8_t *)malloc(N_neurons * N_neurons * sizeof(uint8_t));

        for(alpha = alpha_0; alpha <= alpha_f; alpha += d_alpha){
            
            for (int i = 0; i < N_neurons; i++){
                for (int j = 0; j < N_neurons; j++){
                    AdjMatrix_flat[i*N_neurons + j] = 0;
                    AdjMatrix[i][j] = 0;
                    
                }
            }
            // printf("Patata\n");

            // sprintf(filename, "networks/neurons_params_a%.4f.txt", alpha);
            // FILE *neurons_params;
            // if ((neurons_params = fopen(filename, "w")) == NULL) {
            //     printf("Error opening file %s\n", filename);
            //     exit(1);
            // }
            // fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

            // sprintf(filename, "networks/adjacency_matrix_a%.4f.txt", alpha);

            // FILE *adjacency_matrix;
            // if ((adjacency_matrix = fopen(filename, "w")) == NULL){
            //     printf("Error opening file %s\n", filename);
            //     exit(1);
            // }
            if (aggregation){ // This is in the case of aggregation, we read the positions and parameters from the file created in tune_positions.c
                FILE *neurons_params;
                char filename_np[MAX_STRING_LENGTH];
                sprintf(filename_np, "agrupation/neurons_params_agg1_ng1_nic50_m0.50_s0.10.txt");
                if ((neurons_params = fopen(filename_np, "r")) == NULL) {
                    printf("Error opening file %s\n", filename_np);
                    exit(1);
                }
                char header[MAX_STRING_LENGTH];
                fgets(header, sizeof(header), neurons_params); // Skip header line
                for(int i = 0; i < N_neurons; i++){
                    if (fscanf(neurons_params, "%lf\t%lf\t%lf\t%lf\t%lf\n", &X[i], &Y[i], &soma_diameter, &dendrites_diameters[i], &axon_lengths[i]) != 5) {
                        printf("Error reading line %d from file %s\n", i+2, filename_np);
                        exit(1);
                    }
                }
                fclose(neurons_params);
                printf("Neurons parameters loaded for alpha = %.4f. Creating connections... (Run %d/%d)\n", alpha, run+1, runs);
            }
            else{
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

                    do {
                        dendrites_diameters[i] = d_mean + d_sigma * box_muller();
                    } while (dendrites_diameters[i] <= 0.0);
                    double u = randomInPR(0.0, 1.0);
                    axon_lengths[i] = inverse_cumulative_rayleigh(u, sigma_rayleigh);
                    // fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
                
                }
                printf("Neurons placed for alpha = %.4f. Creating connections... (Run %d/%d)\n", alpha, run+1, runs);
            }
            
    
            
            for(int i = 0; i < N_neurons; i++){
                double segment_vector_x, segment_vector_y;
                double initial_angle = randomInPR(0.0, 2.0*PI);
                segment_vector_x = X[i] + (soma_diameter/2.0) * cos(initial_angle);
                segment_vector_y = Y[i] + (soma_diameter/2.0) * sin(initial_angle);
                int number_of_segments = (int)(axon_lengths[i] / segment_length);
                for(int k = 0; k < number_of_segments; k++){

                    double dx,dy;
                    double angle;
                    double end_segment_x, end_segment_y;

                    angle = initial_angle + box_muller() * sigma_axon_angle; // Standard deviation of sigma_axon_angle radians
                    new_vector_segment(segment_length, angle, &dx, &dy);
                    // printf("Neuron %d, segment %d: angle = %f, dx = %f, dy = %f\n", i, k, angle, dx, dy);
                    end_segment_x = segment_vector_x + dx;
                    end_segment_y = segment_vector_y + dy;
                    PBC(&end_segment_x, L);
                    PBC(&end_segment_y, L);
                    initial_angle = angle;
                    
                    for(int j = 0; j<N_neurons; j++){
                        if (i != j){
                            bool intersection = new_axon_intersection(X, Y, dendrites_diameters, i, j, &segment_vector_x, &segment_vector_y, L, dx, dy);
                            if (intersection){
                                // printf("Intersection detected between neuron %d and neuron %d by segment k = %d\n", i, j, k);
                                number_of_intersections++;
                                if (randomInPR(0.0, 1.0) < alpha){ // Connection probability alpha
                                    AdjMatrix[i][j] = 1;
                                }
                            }
                        }
                    }

                    segment_vector_x = end_segment_x;
                    segment_vector_y = end_segment_y;

                }
                mean_intersections += number_of_intersections;
                number_of_intersections = 0;    
            }

            for(int i = 0; i < N_neurons; i++){
                for(int j = 0; j < N_neurons; j++){
                    // fprintf(adjacency_matrix, "%d", AdjMatrix[i][j]);
                    AdjMatrix_flat[i*N_neurons + j] = (uint8_t)AdjMatrix[i][j];
                }
                // fprintf(adjacency_matrix, "\n");
            }

            fwrite(AdjMatrix_flat, sizeof(uint8_t), N_neurons * N_neurons, fb);  

            // fclose(adjacency_matrix);
            // fclose(neurons_params);
        
        }

        free(AdjMatrix_flat);
        fclose(fb);
    }
    
    printf("Mean intersections per neuron: %f\n", mean_intersections / N_neurons);

    // ? #############################################

    // fclose(test_rayleigh);
    // fclose(test_gaussian); 

    // Al final, libera la memoria
    for (int i = 0; i < NMAX; i++) {
        free(AdjMatrix[i]);
    }
    free(AdjMatrix);
    free(X);
    free(Y);
    free(dendrites_diameters);
    free(axon_lengths);

    return 0;
}