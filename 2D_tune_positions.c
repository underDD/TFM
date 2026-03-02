#include "2D_head.h"

int main(){

    ini_ran(time(NULL));
    
    double *X, *Y; // Neuron positions

    double L, rho; int N_neurons;
    double soma_diameter; 
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh; 
    double segment_length, sigma_axon_angle;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];

    int agg; // 0 no aggragation, 1 aggregation
    
    L = 2.0; // Box length (Everything in mm)
    rho = 300; // Neuron density (neurons/mm^2)
    N_neurons = (int)(rho * L * L);
    
    soma_diameter = 0.015; // Fixed soma diameter (radius of 7.5 microns)
    
    d_mean = 0.15; // Mean of dendrite diameter distribution (Gaussian) (300 microns)
    d_sigma = 0.02; // Standard deviation of dendrite diameter distribution (Gaussian) (40 microns)
 
    l_mean = 1; // Mean of axon length distribution (Rayleigh) (mm)
    sigma_rayleigh = l_mean*sqrt(2/PI); // Mean of axon length distribution (Rayleigh)
    l_sigma = sqrt((4 - PI)/2.0) * sigma_rayleigh; // Calculate standard deviation from sigma

    segment_length = 0.01; // Length of each segment in the axon (10 microns)
    sigma_axon_angle = 0.1; // Standard deviation of axon angle deviation (radians)

    X = (double *)malloc(N_neurons * sizeof(double));
    Y = (double *)malloc(N_neurons * sizeof(double));

    somas = (double *)malloc(N_neurons * sizeof(double));
    dendrites_diameters = (double *)malloc(N_neurons * sizeof(double));
    axon_lengths = (double *)malloc(N_neurons * sizeof(double));

    agg = 0; // Control variable for aggregation (0 no aggregation, 1 aggregation)

    if (agg == 0){

        FILE *neurons_params;
        sprintf(filename, "2D_initial_configurations/2D_neurons_params_random_L%.1lf_rho%.0f_l%.2lf.txt", L, rho, l_mean);
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }
        printf("Placing nurons randomly in the box...\n Prameters:\n");
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

        fprintf(neurons_params, "L = %.1lf\n", L);
        fprintf(neurons_params, "rho = %.0f\n", rho);
        fprintf(neurons_params, "soma_diameter = %.3lf\n", soma_diameter);
        fprintf(neurons_params, "d_mean = %.3lf\n", d_mean);
        fprintf(neurons_params, "d_sigma = %.3lf\n", d_sigma);
        fprintf(neurons_params, "l_mean = %.3lf\n", l_mean);
        fprintf(neurons_params, "l_sigma = %.3lf\n", l_sigma);
        fprintf(neurons_params, "sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        fprintf(neurons_params, "segment_length = %.3lf\n", segment_length);
        fprintf(neurons_params, "sigma_axon_angle = %.3lf\n", sigma_axon_angle);
        fprintf(neurons_params, "nc = %d\n", 0);
        fprintf(neurons_params, "base_sigma = %.3lf\n", 0.0);
        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        bool overlap;
        for(int i = 0; i < N_neurons; i++){ // Only to place the neurons in the plane, without making links yet
            do{

                X[i] = randomInPR(-L/2.0, L/2.0);
                Y[i] = randomInPR(-L/2.0, L/2.0);
                
                overlap = false;
                for(int j = 0; j < i; j++){
                    double dx = min_image(X[i]-X[j], L);
                    double dy = min_image(Y[i]-Y[j], L);
                    double dist = sqrt(dx*dx + dy*dy);
                    if (dist < soma_diameter){
                        overlap = true;
                        break;
                    }
                }

            }while(overlap);

            somas[i] = soma_diameter;

            do {
                dendrites_diameters[i] = d_mean + d_sigma * box_muller();
            } while (dendrites_diameters[i] <= 0.0);

            axon_lengths[i] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
            fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
            printf("Placed neuron %d/%d\r", i+1, N_neurons);
        }

        fclose(neurons_params);

    }

    else {

        int n_centers = (int)(N_neurons * 0.1/4);
        int neurons_in_center = (int)(N_neurons / n_centers);
        double mean_x, mean_y;
        double sigma_x, sigma_y, base_sigma, variation_sigma;

        base_sigma = 0.05; // Base sigma for the Gaussian distribution of the centers of aggregation
        variation_sigma = base_sigma * 0.1; // Variation of sigma for each center (10% of the base sigma)

        sprintf(filename, "2D_initial_configurations/2D_neurons_params_agg_nc%d_s%.2f_L%.1f_rho%.0f_l%.2f.txt", n_centers, base_sigma, L, rho, l_mean);
        FILE *neurons_params;

        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        printf("Placing nurons randomly in the box...\n Prameters:\n");
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

        fprintf(neurons_params, "L = %.1lf\n", L);
        fprintf(neurons_params, "rho = %.0f\n", rho);
        fprintf(neurons_params, "soma_diameter = %.3lf\n", soma_diameter);
        fprintf(neurons_params, "d_mean = %.3lf\n", d_mean);
        fprintf(neurons_params, "d_sigma = %.3lf\n", d_sigma);
        fprintf(neurons_params, "l_mean = %.3lf\n", l_mean);
        fprintf(neurons_params, "l_sigma = %.3lf\n", l_sigma);
        fprintf(neurons_params, "sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        fprintf(neurons_params, "segment_length = %.3lf\n", segment_length);
        fprintf(neurons_params, "sigma_axon_angle = %.3lf\n", sigma_axon_angle);
        fprintf(neurons_params, "nc = %d\n", n_centers);
        fprintf(neurons_params, "base_sigma = %.3lf\n", base_sigma);

        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        for (int c = 0; c < n_centers; c++) {

            sigma_x =  base_sigma + randomInPR(-variation_sigma, variation_sigma); 
            sigma_y =  base_sigma + randomInPR(-variation_sigma, variation_sigma);
            
            mean_x = randomInPR(-L/2.0, L/2.0);
            mean_y = randomInPR(-L/2.0, L/2.0);
            
            for (int n = 0; n < neurons_in_center; n++) {
                bool overlap;
                do {

                    X[c*neurons_in_center + n] = mean_x + sigma_x * box_muller();
                    Y[c*neurons_in_center + n] = mean_y + sigma_y * box_muller();
                    PBC(&X[c*neurons_in_center + n], L);
                    PBC(&Y[c*neurons_in_center + n], L);
                    
                    overlap = false;
                    for (int j = 0; j < c*neurons_in_center + n; j++) {
                        double dx = X[c*neurons_in_center + n] - X[j]; dx = min_image(dx, L);
                        double dy = Y[c*neurons_in_center + n] - Y[j]; dy = min_image(dy, L);
                        double dist = sqrt_distance(dx, dy);
                        if (dist < soma_diameter){ 
                            overlap = true; 
                            break; 
                        }
                    }

                } while (overlap);

                somas[c*neurons_in_center + n] = soma_diameter;
                do { 
                    dendrites_diameters[c*neurons_in_center + n] = d_mean + d_sigma * box_muller(); 
                }while (dendrites_diameters[c*neurons_in_center + n] <= 0.0);

                axon_lengths[c*neurons_in_center + n] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);

                fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[c*neurons_in_center + n], Y[c*neurons_in_center + n], somas[c*neurons_in_center + n], dendrites_diameters[c*neurons_in_center + n], axon_lengths[c*neurons_in_center + n]);
            }   
        } 
        fclose(neurons_params);
    }

    return 0;       
}