#include "3D_head.h"

int main(int argc, char *argv[]) {

    ini_ran(time(NULL));

    if (argc < 7) {
        fprintf(stderr,
            "Usage:\n"
            "%s adjacency_file max_exc_weight max_inh_weight noise_max sim_time seed\n",
            argv[0]
        );
        return 1;
    }

    const char *adj_filename = argv[1];
    double max_exc_weight = atof(argv[2]);
    double max_inh_weight = atof(argv[3]);
    double noise_max = atof(argv[4]);
    int SIM_TIME = atoi(argv[5]);
    unsigned int seed = (unsigned int)atoi(argv[6]);

    srand(seed);

    int N = get_matrix_size(adj_filename);
    if (N <= 0) {
        fprintf(stderr, "Could not determine matrix size.\n");
        return 1;
    }

    uint8_t *A = (uint8_t *)malloc(N * N * sizeof(uint8_t));
    if (!A) {
        fprintf(stderr, "Memory allocation error for adjacency matrix.\n");
        return 1;
    }

    if (load_adjacency_matrix(adj_filename, A, N) != 0) {
        free(A);
        return 1;
    }

    char output_filename[MAX_STRING_LENGTH];
    if (make_output_filename(adj_filename, max_exc_weight, output_filename, sizeof(output_filename)) != 0) {
        free(A);
        return 1;
    }

    if (simulate_izhikevich(
            A,
            N,
            SIM_TIME,
            max_exc_weight,
            max_inh_weight,
            noise_max,
            output_filename
        ) != 0) {
        free(A);
        return 1;
    }

    free(A);
    return 0;
}