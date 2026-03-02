#include "3D_head.h"

int main(){

    ini_ran(time(NULL));

    double *X, *Y, *Z;

    double L, rho; int N;
    double soma_diameter;
    double d_mean, d_sigma;
    double l_mean, l_sigma, sigma_rayleigh;
    double segment_length, sigma_pol, sigma_azi;

    double *somas;
    double *dendrites_diameters;
    double *axon_lengths;

    char filename[MAX_STRING_LENGTH];

    int agg; // 0 no aggragation, 1 aggregation

    L = 2.0;
    rho = 300; // neurons per mm^3
    N = (int)(rho * L * L * L);

    soma_diameter = 0.015; // 15 microns

    d_mean = 0.15; // mean of dendrite diameter distribution (Gaussian) (300 microns)
    d_sigma = 0.02; // standard deviation of dendrite diameter distribution (Gaussian) (20 microns)

    l_mean = 1.; // mean of axon length distribution (Rayleigh) (mm)
    sigma_rayleigh = l_mean*sqrt(2/PI); // Sigma parameter for Rayleigh distribution with given mean
    l_sigma = sqrt((4 - PI)/2.0) * sigma_rayleigh; // Calculate standard deviation from sigma
    
    segment_length = 0.01;
    sigma_pol = 0.1; // Standard deviation of polar angle distribution (Gaussian) (radians)
    sigma_azi = 0.1; // Standard deviation of azimuthal angle distribution (Gaussian) (radians)

    X = (double *)malloc(N * sizeof(double));
    Y = (double *)malloc(N * sizeof(double));
    Z = (double *)malloc(N * sizeof(double));

    somas = (double *)malloc(N * sizeof(double));
    dendrites_diameters = (double *)malloc(N * sizeof(double));
    axon_lengths = (double *)malloc(N * sizeof(double));

    agg = 1; // Control variable for aggregation (0 no aggregation, 1 aggregation)

    if (agg == 0){
        
        FILE *neurons_params;
        sprintf(filename, "3D_initial_configurations/3D_neurons_params_random_L%.1lf_rho%.0f_l%.2lf.txt", L, rho, l_mean);
        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        printf("Placing nurons randomly in the box...\n Prameters:\n");
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

        fprintf(neurons_params, "L = %.1lf\n", L);
        fprintf(neurons_params, "rho = %.0f\n", rho);
        fprintf(neurons_params, "soma_diameter = %.3lf\n", soma_diameter);
        fprintf(neurons_params, "d_mean = %.3lf\n", d_mean);
        fprintf(neurons_params, "d_sigma = %.3lf\n", d_sigma);
        fprintf(neurons_params, "l_mean = %.3lf\n", l_mean);
        fprintf(neurons_params, "l_sigma = %.3lf\n", l_sigma);
        fprintf(neurons_params, "sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        fprintf(neurons_params, "segment_length = %.3lf\n", segment_length);
        fprintf(neurons_params, "sigma_pol = %.3lf\n", sigma_pol);
        fprintf(neurons_params, "sigma_azi = %.3lf\n", sigma_azi);
        fprintf(neurons_params, "nc = %d\n", 0);
        fprintf(neurons_params, "base_sigma = %.3lf\n", 0.0);

        fprintf(neurons_params, "X\tY\tZ\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        bool overlap;
        for(int i = 0; i < N; i++){
            do{

                X[i] = randomInPR(-L/2.0, L/2.0);
                Y[i] = randomInPR(-L/2.0, L/2.0);
                Z[i] = randomInPR(-L/2.0, L/2.0);
                
                overlap = false;
                for(int j = 0; j < i; j++){
                    double dx = min_image(X[i]-X[j], L); // Considerar condiciones de frontera periódicas
                    double dy = min_image(Y[i]-Y[j], L);
                    double dz = min_image(Z[i]-Z[j], L);
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
            printf("Placed neuron %d/%d\r", i+1, N);
        }

        fclose(neurons_params);

    }else {
        int n_centers = (int)(N * 0.1/6); // Number of aggregation centers (10% of total neurons)
        int neurons_in_center = (int)(N / n_centers); // Number of neurons in each center
        double mean_x, mean_y, mean_z, sigma_x, sigma_y, sigma_z;
        double base_sigma, variation_sigma;

        base_sigma = 0.1; // Base sigma for the Gaussian distribution of the centers of aggregation
        variation_sigma = base_sigma * 0.1; // Variation of sigma for each center (10% of the base sigma)

        sprintf(filename, "3D_initial_configurations/3D_neurons_params_agg_nc%d_s%.2lf_L%.1lf_rho%.0f_l%.2lf.txt", n_centers, base_sigma, L, rho, l_mean);
        FILE *neurons_params;

        if ((neurons_params = fopen(filename, "w")) == NULL) {
            printf("Error opening file %s\n", filename);
            exit(1);
        }

        printf("Placing nurons randomly in the box...\n Prameters:\n");
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
        printf("base_sigma = %.3lf\n", base_sigma);

        fprintf(neurons_params, "L = %.1lf\n", L);
        fprintf(neurons_params, "rho = %.0f\n", rho);
        fprintf(neurons_params, "soma_diameter = %.3lf\n", soma_diameter);
        fprintf(neurons_params, "d_mean = %.3lf\n", d_mean);
        fprintf(neurons_params, "d_sigma = %.3lf\n", d_sigma);
        fprintf(neurons_params, "l_mean = %.3lf\n", l_mean);
        fprintf(neurons_params, "l_sigma = %.3lf\n", l_sigma);
        fprintf(neurons_params, "sigma_rayleigh = %.3lf\n", sigma_rayleigh);
        fprintf(neurons_params, "segment_length = %.3lf\n", segment_length);
        fprintf(neurons_params, "sigma_pol = %.3lf\n", sigma_pol);
        fprintf(neurons_params, "sigma_azi = %.3lf\n", sigma_azi);
        fprintf(neurons_params, "nc = %d\n", n_centers);
        fprintf(neurons_params, "s = %.3lf\n", base_sigma);

        fprintf(neurons_params, "X\tY\tZ\tSoma_Diameter\tDendrite_Diameter\tAxon_Length\n");

        for(int c = 0; c < n_centers; c++){

            sigma_x = base_sigma + randomInPR(-variation_sigma, variation_sigma);
            sigma_y = base_sigma + randomInPR(-variation_sigma, variation_sigma);
            sigma_z = base_sigma + randomInPR(-variation_sigma, variation_sigma);

            mean_x = randomInPR(-L/2.0, L/2.0);
            mean_y = randomInPR(-L/2.0, L/2.0);
            mean_z = randomInPR(-L/2.0, L/2.0);

            for(int n = 0; n < neurons_in_center; n++){
                bool overlap;
                do{

                    X[c*neurons_in_center + n] = mean_x + sigma_x * box_muller();
                    Y[c*neurons_in_center + n] = mean_y + sigma_y * box_muller();
                    Z[c*neurons_in_center + n] = mean_z + sigma_z * box_muller();
                    PBC(&X[c*neurons_in_center + n], L);
                    PBC(&Y[c*neurons_in_center + n], L);
                    PBC(&Z[c*neurons_in_center + n], L);

                    overlap = false;
                    for(int j = 0; j < c*neurons_in_center + n; j++){
                        double dx = min_image(X[c*neurons_in_center + n]-X[j], L);
                        double dy = min_image(Y[c*neurons_in_center + n]-Y[j], L);
                        double dz = min_image(Z[c*neurons_in_center + n]-Z[j], L);
                        double dist = sqrt(dx*dx + dy*dy + dz*dz);
                        if (dist < soma_diameter){
                            overlap = true;
                            break;
                        }
                    }

                }while(overlap);

                somas[c*neurons_in_center + n] = soma_diameter;
                do{
                    dendrites_diameters[c*neurons_in_center + n] = d_mean + d_sigma * box_muller();
                }while(dendrites_diameters[c*neurons_in_center + n] <= 0.0);

                axon_lengths[c*neurons_in_center + n] = inverse_cumulative_rayleigh(randomInPR(0.0, 1.0), sigma_rayleigh);
                fprintf(neurons_params, "%f\t%f\t%f\t%f\t%f\t%f\n", X[c*neurons_in_center + n], Y[c*neurons_in_center + n], Z[c*neurons_in_center + n], somas[c*neurons_in_center + n], dendrites_diameters[c*neurons_in_center + n], axon_lengths[c*neurons_in_center + n]);
            }
        }
        fclose(neurons_params);
    }

    return 0;
}
