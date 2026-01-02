#include "head.h"

int main(){
        printf("patata\n");

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

    double ag;
    ag = 0.0; // Level of aggrupation 

    printf("Dendrite diameter distribution (Gaussian): mean = %f, sigma = %f\n", d_mean, d_sigma);
    printf("Axon length distribution (Rayleigh): mean = %f, sigma = %f\n", l_mean, l_sigma);

    X = (double *)malloc(N_neurons * sizeof(double));
    Y = (double *)malloc(N_neurons * sizeof(double));
    dendrites_diameters = (double *)malloc(N_neurons * sizeof(double));
    axon_lengths = (double *)malloc(N_neurons * sizeof(double));

    sprintf(filename, "agrupation/neurons_params_agg%.4f.txt", ag);
    FILE *neurons_params;
    if ((neurons_params = fopen(filename, "w")) == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

    for(int i = 0; i < N_neurons; i++){
        X[i] = Y[i] = 0.0;
        dendrites_diameters[i] = 0.0;
        axon_lengths[i] = 0.0;
    }

    int gaussian_clusters = 4;
    int neurons_in_cluster = 50;
    double mean_x[gaussian_clusters] = {-0.5, 0.5, 0.5, -0.5};
    double mean_y[gaussian_clusters] = {0.5, 0.5, -0.5, -0.5};
    double sigma_x[gaussian_clusters] = {0.1, 0.1, 0.1, 0.1};
    double sigma_y[gaussian_clusters] = {0.1, 0.1, 0.1, 0.1};


    

    if (ag == 0.0){
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
    }

    else {
        bool overlap;
        for(int i = 0; i < N_neurons; i++){ // Only to place the neurons in the plane, without making links yet
            do{
                
                overlap = false;
                if (i < 50){
                    X[i] = mean_x[0] + sigma_x[0] * box_muller();
                    Y[i] = mean_y[0] + sigma_y[0] * box_muller();
                }
                else if (i < 100 && i >= 50){
                    X[i] = mean_x[1] + sigma_x[1] * box_muller();
                    Y[i] = mean_y[1] + sigma_y[1] * box_muller();
                }
                else if (i < 150 && i >= 100){
                    X[i] = mean_x[2] + sigma_x[2] * box_muller();
                    Y[i] = mean_y[2] + sigma_y[2] * box_muller();
                }
                else if (i < 200 && i >= 150){
                    X[i] = mean_x[3] + sigma_x[3] * box_muller();
                    Y[i] = mean_y[3] + sigma_y[3] * box_muller();
                }
                else {
                    X[i] = randomInPR(-L/2.0, L/2.0);
                    Y[i] = randomInPR(-L/2.0, L/2.0);
                }

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
    }
    
    return 0;
}