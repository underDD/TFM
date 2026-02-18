#include "head.h"

int main(){
        printf("patata\n");

    // ? ############### VARIABLES ###############

    double l_mean, l_sigma; 
    double d_mean, d_sigma;

    int N_neurons;
    double soma_diameter; 

    double L; 

    double *X, *Y; // Neuron positions
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];

    
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
    l_mean = 0.5; // Mean of axon length distribution (Rayleigh) (mm)
    // l_sigma = sqrt(0.8); // Sigma parameter for Rayleigh distribution with given mean
    sigma_rayleigh = l_mean*sqrt(2/PI); // Mean of axon length distribution (Rayleigh)
    // sigma_rayleigh = l_sigma*sqrt(2/(4 - PI)); // Recalculate sigma from desired standard deviation
    // l_mean = sigma_rayleigh*sqrt(PI/2.0); // Calculate mean from sigma
    l_sigma = sqrt((4 - PI)/2.0) * sigma_rayleigh; // Calculate standard deviation from sigma

    int agg = 1; // 0 no aggragation, 1 aggregation

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

    if (agg == 0){

        sprintf(filename, "configurations/neurons_params_random_l%.2lf.txt", l_mean);
        FILE *neurons_params;
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }
        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

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
            fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n", X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
        }
        fclose(neurons_params);

    }

    else {

        int n_centers = 25;
        int neurons_in_center = (int)(N_neurons / n_centers);
        double mean_x, mean_y;
        double sigma_x, sigma_y, base_sigma;

        base_sigma = 0.25; // Base sigma for the Gaussian distribution of the centers of aggregation

        sprintf(filename, "configurations/neurons_params_agg_nc%d_s%.2f.txt", n_centers, base_sigma);
        FILE *neurons_params;

        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        for (int c = 0; c < n_centers; c++) {
            // Here I define the parameters of the Gaussian distribution of the centar of aggregation. sigma = 0.1 + little variation
            sigma_x = sigma_y = base_sigma + randomInPR(-0.02, 0.02); 
            mean_x = randomInPR(-L/2.0, L/2.0);
            mean_y = randomInPR(-L/2.0, L/2.0);
            for (int n = 0; n < neurons_in_center; n++) {

                bool overlap;

                do {

                    overlap = false;

                    X[c*neurons_in_center + n] = mean_x + sigma_x * box_muller();
                    Y[c*neurons_in_center + n] = mean_y + sigma_y * box_muller();
                    PBC(&X[c*neurons_in_center + n], L);
                    PBC(&Y[c*neurons_in_center + n], L);
                    
                    for (int j = 0; j < c*neurons_in_center + n; j++) {

                        double dx = X[c*neurons_in_center + n] - X[j]; dx = min_image(dx, L);
                        double dy = Y[c*neurons_in_center + n] - Y[j]; dy = min_image(dy, L);
                        double dist = sqrt_distance(dx, dy);

                        if (dist < soma_diameter) { overlap = true; break; }
                    }

                } while (overlap);

                // resto de parámetros
                do { 
                    dendrites_diameters[c*neurons_in_center + n] = d_mean + d_sigma * box_muller(); 
                }
                while (dendrites_diameters[c*neurons_in_center + n] <= 0.0);

                double u = randomInPR(0.0, 1.0);
                axon_lengths[c*neurons_in_center + n] = inverse_cumulative_rayleigh(u, sigma_rayleigh);

                fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n",
                        X[c*neurons_in_center + n], Y[c*neurons_in_center + n], soma_diameter, dendrites_diameters[c*neurons_in_center + n], axon_lengths[c*neurons_in_center + n]);
            }   
        } 

        fclose(neurons_params);
    }

    return 0;       
}