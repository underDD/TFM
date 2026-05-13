#include "2D_head.h"

int main(int argc, char *argv[]) {

    ini_ran(time(NULL));
    
    double *X, *Y; // Neuron positions

    double L, rho; int N_neurons;
    double soma_diameter; 
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh; 
    double segment_length, sigma_axon_angle;
    double R;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];
    char BC_type[MAX_STRING_LENGTH];
    char geometry[MAX_STRING_LENGTH];

    int agg; // 0 no aggragation, 1 aggregation
    int n_centers;
    double base_sigma;

    if (argc != 14) {
        printf("Uso:\n%s L rho soma_diameter d_mean d_sigma l_mean segment_length sigma_axon_angle agg n_centers base_sigma PBC geometry\n", argv[0]);
        return 1;
    }

    L = atof(argv[1]); // Could be the length of the square or the diameter of the circle
    rho = atof(argv[2]);
    soma_diameter = atof(argv[3]);
    d_mean = atof(argv[4]);
    d_sigma = atof(argv[5]);
    l_mean = atof(argv[6]);
    segment_length = atof(argv[7]);
    sigma_axon_angle = atof(argv[8]);
    agg = atoi(argv[9]);
    n_centers = atoi(argv[10]);
    base_sigma = atof(argv[11]);
    strcpy(BC_type, argv[12]);
    strcpy(geometry, argv[13]);

    if (strcmp(geometry, "square") == 0){
        N_neurons = (int)(rho * L * L);
    }
    else if (strcmp(geometry, "circle") == 0){
        R = L / 2; // Calculate the radius of the circle to maintain the same area as the square
        N_neurons = (int)(rho * PI * R * R);
    }
    
    sigma_rayleigh = l_mean*sqrt(2/PI); // Mean of axon length distribution (Rayleigh)
    l_sigma = sqrt((4 - PI)/2.0) * sigma_rayleigh; // Calculate standard deviation from sigma

    // strcpy(BC_type, "PBC"); // Boundary conditions type (PBC: Periodic Boundary Conditions, RBC: Reflective Boundary Conditions)

    X = (double *)malloc(N_neurons * sizeof(double));
    Y = (double *)malloc(N_neurons * sizeof(double));

    somas = (double *)malloc(N_neurons * sizeof(double));
    dendrites_diameters = (double *)malloc(N_neurons * sizeof(double));
    axon_lengths = (double *)malloc(N_neurons * sizeof(double));

    
    if (agg == 0){
        
        FILE *neurons_params;
        sprintf(filename, "2D_initial_configurations/2D_neurons_params_%s_%s_random_L%.1lf_rho%.0f_d%.2f.txt", geometry, BC_type, L, rho, d_mean);
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }
        
        printf("ESTOY EJECUTANDO LA VERSION NUEVA\n");
        printf("Placing neurons randomly...\nParameters:\n");
        printf("L = %.1lf\n", L);
        printf("rho = %.0f\n", rho);
        printf("N = %d\n", N_neurons);
        printf("soma_diameter = %.3lf\n", soma_diameter);
        printf("d_mean = %.3lf\n", d_mean);
        printf("d_sigma = %.3lf\n", d_sigma);
        printf("l_mean = %.3lf\n", l_mean);
        printf("l_sigma = %.3lf\n", l_sigma);
        printf("sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        printf("segment_length = %.3lf\n", segment_length);
        printf("sigma_axon_angle = %.3lf\n", sigma_axon_angle);
        printf("Geometry = %s\n", geometry);
        printf("Output file: %s\n", filename);  
        
        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");
        
        bool overlap;
        for(int i = 0; i < N_neurons; i++){ // Only to place the neurons in the plane, without making links yet
            do{
                
                if (strcmp(geometry, "square") == 0){
                    double margin = soma_diameter/2.0; // Margin to avoid placing neurons too close to the borders when using non-periodic BC
                    X[i] = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    Y[i] = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                }
                else if (strcmp(geometry, "circle") == 0){
                    double R_eff = R - soma_diameter/2.0; // Effective radius to ensure the entire soma fits within the circle
                    double radius = R_eff * sqrt(randomInPR(0.0, 1.0)); // Random radius with uniform distribution in the circle
                    double angle = randomInPR(0.0, 2*PI); // Random angle
                    X[i] = radius * cos(angle); 
                    Y[i] = radius * sin(angle);
                }
                
                overlap = false;
                double dx = 0, dy = 0;
                for(int j = 0; j < i; j++){

                    if (strcmp(BC_type, "PBC") == 0){
                        dx = min_image(X[i]-X[j], L);
                        dy = min_image(Y[i]-Y[j], L);
                    }
                    else{
                        dx = X[i]-X[j];
                        dy = Y[i]-Y[j];
                    }

                    double dist = sqrt(dx*dx + dy*dy);
                    if (dist < soma_diameter){
                        overlap = true;
                        break;
                    }
                }
                
            }while(overlap);
            
            somas[i] = soma_diameter;
            
            do{
                dendrites_diameters[i] = d_mean + d_sigma * box_muller();
            } while (dendrites_diameters[i] <= 0.0);
            
            axon_lengths[i] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
            fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
            // printf("Placed neuron %d/%d\r", i+1, N_neurons);
            fflush(stdout);
        }
        printf("\n");
        
        fclose(neurons_params);
        
    }
    
    else {
        
        int base_count = N_neurons / n_centers;
        int remainder = N_neurons % n_centers;
        double mean_x, mean_y;
        double sigma_x, sigma_y, variation_sigma;
        
        // base_sigma = 0.05; // Base sigma for the Gaussian distribution of the centers of aggregation
        variation_sigma = base_sigma * 0.1; // Variation of sigma for each center (10% of the base sigma)
        
        sprintf(filename, "2D_initial_configurations/2D_neurons_params_%s_%s_agg_nc%d_s%.4f_L%.1f_rho%.0f_d%.2f.txt", geometry, BC_type, n_centers, base_sigma, L, rho, d_mean);
        FILE *neurons_params;
        
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }
        
        printf("Placing neurons with aggregation...\n Parameters:\n");
        printf("L = %.1lf\n", L);
        printf("rho = %.0f\n", rho);
        printf("N = %d\n", N_neurons);
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
        printf("Geometry = %s\n", geometry);
        printf("Output file: %s\n", filename);

        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");
        
        int idx = 0;

        for (int c = 0; c < n_centers; c++) {

            int n_this_center = base_count + (c < remainder ? 1 : 0); // Distribute the remainder neurons among the first 'remainder' centers
            
            sigma_x =  base_sigma + randomInPR(-variation_sigma, variation_sigma); 
            sigma_y =  base_sigma + randomInPR(-variation_sigma, variation_sigma);
            
            if (strcmp(BC_type, "PBC") == 0){
                mean_x = randomInPR(-L/2.0, L/2.0);
                mean_y = randomInPR(-L/2.0, L/2.0);
            }
            else{

                if (strcmp(geometry, "square") == 0){
                    double margin = soma_diameter/2.0; // Margin to avoid placing centers too close to the borders when using non-periodic BC            
                    mean_x = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                    mean_y = randomInPR(-L/2.0 + margin, L/2.0 - margin);
                }
                else if (strcmp(geometry, "circle") == 0){
                    double margin = soma_diameter/2.0; // Margin to avoid placing centers too close to the borders when using non-periodic BC
                    double R_eff = R - margin; // Effective radius to ensure centers are not too close to the border
                    double radius = R_eff * sqrt(randomInPR(0.0, 1.0)); // Random radius with uniform distribution in the circle
                    double angle = randomInPR(0.0, 2*PI); // Random angle
                    mean_x = radius * cos(angle); 
                    mean_y = radius * sin(angle);
                }
                
            }
            
            for (int n = 0; n < n_this_center; n++) {
                bool overlap;
                do {
                    
                    if (strcmp(BC_type, "PBC") == 0){
                        X[idx] = mean_x + sigma_x * box_muller();
                        Y[idx] = mean_y + sigma_y * box_muller();
                        PBC(&X[idx], L);
                        PBC(&Y[idx], L);
                    }
                    else{
                        double margin = soma_diameter / 2.0; // Margin to avoid placing neurons too close to the walls
                        if (strcmp(geometry, "square") == 0){
                            do {
                                X[idx] = mean_x + sigma_x * box_muller();
                                Y[idx] = mean_y + sigma_y * box_muller();
                            } while (X[idx] < -L/2.0 + margin || X[idx] > L/2.0 - margin ||
                                    Y[idx] < -L/2.0 + margin || Y[idx] > L/2.0 - margin);
                        }
                        else if (strcmp(geometry, "circle") == 0){
                            double margin = soma_diameter/2.0; // Margin to avoid placing centers too close to the borders when using non-periodic BC
                            double R_eff = R - margin;
                            do{
                                X[idx] = mean_x + sigma_x * box_muller();
                                Y[idx] = mean_y + sigma_y * box_muller();
                            }while(X[idx]*X[idx] + Y[idx]*Y[idx] > R_eff*R_eff); // Ensure the point is within the circle
                        }
                    }
                        

                    overlap = false;
                    double dx = 0, dy = 0;

                    for (int j = 0; j < idx; j++) {
                        if (strcmp(BC_type, "PBC") == 0){
                            dx = min_image(X[idx] - X[j], L);
                            dy = min_image(Y[idx] - Y[j], L);
                        }
                        else{
                            dx = X[idx] - X[j];
                            dy = Y[idx] - Y[j];
                        }
                        double dist = sqrt_distance(dx, dy);
                        if (dist < soma_diameter){ 
                            overlap = true; 
                            break; 
                        }
                    }
                    
                } while (overlap);
                
                somas[idx] = soma_diameter;
                do { 
                    dendrites_diameters[idx] = d_mean + d_sigma * box_muller(); 
                }while (dendrites_diameters[idx] <= 0.0);
                
                axon_lengths[idx] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
                
                fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[idx], Y[idx], somas[idx], dendrites_diameters[idx], axon_lengths[idx]);
                // printf("Placed neuron %d/%d\r", idx + 1, N_neurons);
                fflush(stdout);
                idx++;
            }   
            printf("\n");
        } 
        fclose(neurons_params);
    }
    
    free(X);
    free(Y);
    free(somas);
    free(dendrites_diameters);
    free(axon_lengths);
    return 0;       
}



// rho = 250; // Neuron density (neurons/mm^2)
// L = 3.0; // Box length (Everything in mm)
// soma_diameter = 0.015; // Fixed soma diameter (radius of 7.5 microns)
// d_mean = 0.2; // Mean of dendrite diameter distribution (Gaussian) (300 microns)
// d_sigma = 0.02; // Standard deviation of dendrite diameter distribution (Gaussian) (40 microns)
// l_mean = 0.5; // Mean of axon length distribution (Rayleigh) (mm)
// segment_length = 0.01; // Length of each segment in the axon (10 microns)
// sigma_axon_angle = 0.1; // Standard deviation of axon angle deviation (radians)
// agg = 1; // Control variable for aggregation (0 no aggregation, 1 aggregation)

// if (strcmp(geometry, "circle") == 0){
//         bool overlap;
//         double angle, radius;
//         for(int i = 0; i < N_neurons; i++){
//             do{
//                 radius = randomInPR(0.0, r);
//                 angle = randomInPR(0.0, 2*PI);
//                 overlap = false;

//                 X[i] = radius * cos(angle); 
//                 Y[i] = radius * sin(angle);

//                 for(int j = 0; j < i; j++){
//                     double dx = min_image(X[i]-X[j], L);
//                     double dy = min_image(Y[i]-Y[j], L);
//                     double dist = sqrt(dx*dx + dy*dy);
//                     if (dist < soma_diameter){
//                         overlap = true;
//                         break;
//                     }
//                 }
//             }while(overlap);
//         }
//     }