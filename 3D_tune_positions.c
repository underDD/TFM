#include "3D_head.h"

int main(int argc, char *argv[]) {

    ini_ran(time(NULL));

    double *X, *Y, *Z;

    double L, rho; int N;
    double soma_diameter;
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh;
    double segment_length, sigma_pol, sigma_azi;
    double R, height;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];
    char BC_type[MAX_STRING_LENGTH];
    char geometry[MAX_STRING_LENGTH];

    int agg; // 0 no aggragation, 1 aggregation
    int n_centers;
    double base_sigma;

    if (argc != 15) {
        printf("Uso:\n");
        printf("%s L rho soma_diameter d_mean d_sigma l_mean segment_length sigma_pol sigma_azi agg n_centers base_sigma geometry\n", argv[0]);
        return 1;
    }
 
    L = atof(argv[1]);
    rho = atof(argv[2]);
    soma_diameter = atof(argv[3]);
    d_mean = atof(argv[4]);
    d_sigma = atof(argv[5]);
    l_mean = atof(argv[6]);
    segment_length = atof(argv[7]);
    sigma_pol = atof(argv[8]);
    sigma_azi = atof(argv[9]);
    agg = atoi(argv[10]);
    n_centers = atoi(argv[11]);
    base_sigma = atof(argv[12]);
    strcpy(BC_type, argv[13]);
    strcpy(geometry, argv[14]);
    

    if (strcmp(geometry, "cube") == 0){
        N = (int)(rho * L * L * L);
    }
    else if (strcmp(geometry, "cylinder") == 0){
        R = L / sqrt(PI); // Calculate the radius of the circle to maintain the same area as the square
        height = L; // Keep the same height as the side of the square
        N = (int)(rho * PI * R * R * height);
    }

    // N = (int)(rho * L * L * L);
    sigma_rayleigh = l_mean*sqrt(2/PI); // Sigma parameter for Rayleigh distribution with given mean
    l_sigma = sqrt((4 - PI)/2.0) * sigma_rayleigh; // Calculate standard deviation from sigma
    
    // strcpy(BC_type, "PBC"); // Boundary conditions type (PBC: Periodic Boundary Conditions, RBC: Reflective Boundary Conditions)

    X = (double *)malloc(N * sizeof(double));
    Y = (double *)malloc(N * sizeof(double));
    Z = (double *)malloc(N * sizeof(double));

    somas = (double *)malloc(N * sizeof(double));
    dendrites_diameters = (double *)malloc(N * sizeof(double));
    axon_lengths = (double *)malloc(N * sizeof(double));



    if (agg == 0){
        
        FILE *neurons_params;
        sprintf(filename, "3D_initial_configurations/3D_neurons_params_%s_%s_random_L%.1lf_rho%.0f_d%.2f.txt", geometry, BC_type, L, rho, d_mean);
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        printf("Placing neurons randomly in the box...\nParameters:\n");
        printf("L = %.1lf\n", L);
        printf("rho = %.0f\n", rho);
        printf("N = %d\n", N);
        printf("soma_diameter = %.3lf\n", soma_diameter);
        printf("d_mean = %.3lf\n", d_mean);
        printf("d_sigma = %.3lf\n", d_sigma);
        printf("l_mean = %.3lf\n", l_mean);
        printf("l_sigma = %.3lf\n", l_sigma);
        printf("sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        printf("segment_length = %.3lf\n", segment_length);
        printf("sigma_pol = %.3lf\n", sigma_pol);
        printf("sigma_azi = %.3lf\n", sigma_azi);

        fprintf(neurons_params, "X\tY\tZ\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        double margin = soma_diameter / 2.0; // Margin to avoid placing neurons too close to the walls
        bool overlap;
        for(int i = 0; i < N; i++){
            do{
                if (strcmp(geometry, "cube") == 0){
                    X[i] = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    Y[i] = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    Z[i] = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                }
                else if (strcmp(geometry, "cylinder") == 0){
                    double R_eff = R - margin; // Effective radius to ensure centers are not too close to the border
                    double radius = R_eff * sqrt(randomInPR(0.0, 1.0)); 
                    double angle = randomInPR(0.0, 2*PI); // Random
                    X[i] = radius * cos(angle);
                    Y[i] = radius * sin(angle);
                    Z[i] = randomInPR(-height/2.0 + margin, height/2.0 - margin);
                }
                

                overlap = false;
                double dx = 0, dy = 0, dz = 0;
                for(int j = 0; j < i; j++){

                    if (strcmp(BC_type, "PBC") == 0){
                        dx = min_image(X[i]-X[j], L);
                        dy = min_image(Y[i]-Y[j], L);
                        dz = min_image(Z[i]-Z[j], L);
                    }
                    else{
                        dx = X[i] - X[j];
                        dy = Y[i] - Y[j];
                        dz = Z[i] - Z[j];
                    }

                    double dist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dist < soma_diameter){
                        overlap = true;
                        break;
                    }
                }

            }while(overlap);

            somas[i] = soma_diameter;

            do{
                dendrites_diameters[i] = d_mean + d_sigma * box_muller();
            }while(dendrites_diameters[i] <= 0.0);

            axon_lengths[i] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
            fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], Z[i], somas[i], dendrites_diameters[i], axon_lengths[i]);
            // printf("Placed neuron %d/%d\r", i+1, N);
            fflush(stdout);
        }
        printf("\n");
        fclose(neurons_params);
    }
    
    else {

        int base_count = N / n_centers;
        int remainder = N % n_centers;
        double mean_x, mean_y, mean_z, sigma_x, sigma_y, sigma_z;
        double variation_sigma;

        // base_sigma = 0.1; // Base sigma for the Gaussian distribution of the centers of aggregation
        variation_sigma = base_sigma * 0.1; // Variation of sigma for each center (10% of the base sigma)

        sprintf(filename, "3D_initial_configurations/3D_neurons_params_%s_%s_agg_nc%d_s%.4lf_L%.1lf_rho%.0f_d%.2f.txt", geometry, BC_type, n_centers, base_sigma, L, rho, d_mean);
        FILE *neurons_params;

        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        printf("Placing neurons randomly in the box...\nParameters:\n");
        printf("L = %.1lf\n", L);
        printf("rho = %.0f\n", rho);
        printf("N = %d\n", N);
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
        printf("Output file: %s", filename);

        fprintf(neurons_params, "X\tY\tZ\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        int idx = 0;

        for(int c = 0; c < n_centers; c++){

            int n_this_center = base_count + (c < remainder ? 1 : 0); // Distribute the remainder neurons among the first 'remainder' centers

            sigma_x = base_sigma + randomInPR(-variation_sigma, variation_sigma);
            sigma_y = base_sigma + randomInPR(-variation_sigma, variation_sigma);
            sigma_z = base_sigma + randomInPR(-variation_sigma, variation_sigma);

            if (strcmp(BC_type, "PBC") == 0){
                mean_x = randomInPR(-L/2.0, L/2.0);
                mean_y = randomInPR(-L/2.0, L/2.0);
                mean_z = randomInPR(-L/2.0, L/2.0);
            }
            else{

                if (strcmp(geometry, "cube") == 0){

                    double margin = soma_diameter / 2.0; // Margin to avoid placing neurons too close to the walls
                    mean_x = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    mean_y = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    mean_z = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                }

                else if (strcmp(geometry, "cylinder") == 0){
                    double margin = soma_diameter / 2.0; // Margin to avoid placing neurons too close to the walls
                    double R_eff = R - margin; // Effective radius to ensure centers are not too close to the border
                    double radius = R_eff * sqrt(randomInPR(0.0, 1.0)); 
                    double angle = randomInPR(0.0, 2*PI); // Random
                    mean_x = radius * cos(angle);
                    mean_y = radius * sin(angle);
                    mean_z = randomInPR(-height/2.0 + margin, height/2.0 - margin);
                }
                
            }

            for(int n = 0; n < n_this_center; n++){
                bool overlap;
                do{

                    if (strcmp(BC_type, "PBC") == 0){
                        X[idx] = mean_x + sigma_x * box_muller();
                        Y[idx] = mean_y + sigma_y * box_muller();
                        Z[idx] = mean_z + sigma_z * box_muller();
                        PBC(&X[idx], L);
                        PBC(&Y[idx], L);
                        PBC(&Z[idx], L);
                    }
                    else{

                        if (strcmp(geometry, "cube") == 0){
                            double margin = soma_diameter / 2.0; // Margin to avoid placing neurons too close to the walls
                            do{
                                X[idx] = mean_x + sigma_x * box_muller();
                                Y[idx] = mean_y + sigma_y * box_muller();
                                Z[idx] = mean_z + sigma_z * box_muller();
                            }while(X[idx] < -L/2.0 + margin || X[idx] > L/2.0 - margin ||
                                   Y[idx] < -L/2.0 + margin || Y[idx] > L/2.0 - margin ||
                                   Z[idx] < -L/2.0 + margin || Z[idx] > L/2.0 - margin);
                        }
                        else if (strcmp(geometry, "cylinder") == 0){
                            double margin = soma_diameter / 2.0;
                            double R_eff = R - margin;
                            double radius;
                            double angle;
                            do{
                                X[idx] = mean_x + sigma_x * box_muller();
                                Y[idx] = mean_y + sigma_y * box_muller();
                                Z[idx] = mean_z + sigma_z * box_muller();
                                radius = sqrt(X[idx]*X[idx] + Y[idx]*Y[idx]);
                            }while(radius > R_eff || Z[idx] < -height/2.0 + margin || Z[idx] > height/2.0 - margin);
                        }
                        
                    }

                    overlap = false;
                    double dx = 0, dy = 0, dz = 0;

                    for(int j = 0; j < idx; j++){

                        if (strcmp(BC_type, "PBC") == 0){
                            dx = min_image(X[idx]-X[j], L);
                            dy = min_image(Y[idx]-Y[j], L);
                            dz = min_image(Z[idx]-Z[j], L);
                        }
                        else{
                            dx = X[idx] - X[j];
                            dy = Y[idx] - Y[j];
                            dz = Z[idx] - Z[j];
                        }
                        double dist = sqrt(dx*dx + dy*dy + dz*dz);
                        if (dist < soma_diameter){
                            overlap = true;
                            break;
                        }
                    }

                }while(overlap);

                somas[idx] = soma_diameter;
                do{
                    dendrites_diameters[idx] = d_mean + d_sigma * box_muller();
                }while(dendrites_diameters[idx] <= 0.0);

                axon_lengths[idx] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
                fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\t%f\n", X[idx], Y[idx], Z[idx], somas[idx], dendrites_diameters[idx], axon_lengths[idx]);
                // printf("Placed neuron %d/%d\r", idx + 1, N);
                fflush(stdout);
            }
            printf("\n");
        }
        fclose(neurons_params);

    }

    return 0;
}

// L = 2.0;
// rho = 125; // neurons per mm^3
// soma_diameter = 0.015; // 15 microns
// d_mean = 0.3; // mean of dendrite diameter distribution (Gaussian) (300 microns)
// d_sigma = 0.04; // standard deviation of dendrite diameter distribution (Gaussian) (20 microns)
// l_mean = 1.; // mean of axon length distribution (Rayleigh) (mm)
// segment_length = 0.01;
// sigma_pol = 0.1; // Standard deviation of polar angle distribution (Gaussian) (radians)
// sigma_azi = 0.1; // Standard deviation of azimuthal angle distribution (Gaussian) (radians)
// agg = 1; // Control variable for aggregation (0 no aggregation, 1 aggregation)
