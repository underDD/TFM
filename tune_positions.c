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
    l_mean = 1.0;
    l_sigma = sqrt(0.8); // Sigma parameter for Rayleigh distribution with given mean
    // sigma_rayleigh = l_mean*sqrt(2/PI); // Mean of axon length distribution (Rayleigh)
    sigma_rayleigh = l_sigma*sqrt(2/(4 - PI)); // Recalculate sigma from desired standard deviation
    l_mean = sigma_rayleigh*sqrt(PI/2.0); // Calculate mean from sigma

    int ag = 1; // 0 no aggragation, 1 aggregation

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

    if (ag == 0){

        sprintf(filename, "agrupation/neurons_params_agg%d.txt", ag);
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

        int gaussian_clusters = 4;
        int neurons_in_cluster = 50;
        double mean_x[] = {-0.5, 0.5, 0.5, -0.5};
        double mean_y[] = {0.5, 0.5, -0.5, -0.5};
        double sigma_x[] = {0.1, 0.1, 0.1, 0.1};
        double sigma_y[] = {0.1, 0.1, 0.1, 0.1};

        sprintf(filename, "agrupation/neurons_params_agg%d_ngauss%d_nic%d.txt", ag, gaussian_clusters, neurons_in_cluster);
        FILE *neurons_params;
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }
        fprintf(neurons_params, "X\tY\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        int n_clustered = gaussian_clusters * neurons_in_cluster;

        for (int i = 0; i < N_neurons; i++) {

            int c = -1;
            if (i < n_clustered) c = i / neurons_in_cluster;  // 0..gaussian_clusters-1

            bool overlap;
            do {
                overlap = false;

                if (c >= 0) {
                    X[i] = mean_x[c] + sigma_x[c] * box_muller();
                    Y[i] = mean_y[c] + sigma_y[c] * box_muller();
                    PBC(&X[i], L);
                    PBC(&Y[i], L);
                } else {
                    X[i] = randomInPR(-L/2.0, L/2.0);
                    Y[i] = randomInPR(-L/2.0, L/2.0);
                }

                for (int j = 0; j < i; j++) {
                    // si usas PBC, usa mínima imagen aquí también
                    double dx = X[i] - X[j]; dx = min_image(dx, L);
                    double dy = Y[i] - Y[j]; dy = min_image(dy, L);
                    double dist = sqrt_distance(dx, dy);

                    if (dist < soma_diameter) { overlap = true; break; }
                }

            } while (overlap);

            // resto de parámetros
            do { dendrites_diameters[i] = d_mean + d_sigma * box_muller(); }
            while (dendrites_diameters[i] <= 0.0);

            double u = randomInPR(0.0, 1.0);
            axon_lengths[i] = inverse_cumulative_rayleigh(u, sigma_rayleigh);

            fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\n",
                    X[i], Y[i], soma_diameter, dendrites_diameters[i], axon_lengths[i]);
        }
    }

    return 0;       
}