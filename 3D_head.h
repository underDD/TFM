#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #define MKDIR(path) mkdir(path, 0777)
#endif

#define MAX_STRING_LENGTH 100000
#define fran rand()/((double)RAND_MAX+1)
#define PI 3.14159265358979323846
#define EPS 1E-5
#define NMAX 10000

//**********PARISI RAPUANO*************
#define NormRANu (2.3283063671E-10F)

extern unsigned int irr[256];
extern unsigned int ir1;
extern unsigned char ind_ran,ig1,ig2,ig3;

//************************************

//**********3D_functions.c*************  

float randomIn(double min, double max);
float Random(void);
void ini_ran(int SEMILLA);
float randomInPR(double min, double max);
double box_muller(void);
void PBC(double *x, double L);
double min_image(double d, double L);
void new_vector_segment_3D(double segment_length, double angle_azi, double angle_pol, double *dx, double *dy, double *dz);
bool new_axon_intersection_3D(double *x, double *y, double *z, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double *segment_vector_z, double L, double dx, double dy, double dz);
double sqrt_distance(double x, double y);
double scalar_product(double x1, double y1, double x2, double y2);
double inverse_cumulative_rayleigh(double p, double sigma);
bool new_axon_intersection_noPBC_3D(double cx, double cy, double cz, double dendrite_diameter,
                                double x0, double y0, double z0,
                                double dx, double dy, double dz);
void sticky_walls3D(double initial_segment_vector_x, double initial_segment_vector_y, double initial_segment_vector_z, 
                    double *dx, double *dy, double *dz, double L,
                    double *end_segment_vector_x, double *end_segment_vector_y, double *end_segment_vector_z,
                    double *X, double *Y, double *Z, double *dendrites_diameters,
                    int i, int N_neurons, uint8_t *AdjMatrix_flat,
                    int *trials, int *links, FILE *axon_simulation, int neuron_idx, int segment_idx);

int get_matrix_size(const char *filename);
int load_adjacency_matrix(const char *filename, uint8_t *AdjMatrix_flat, int N);

void init_izhikevich_parameters(
    int N, int Ne, int Ni,
    double *a, double *b, double *c, double *d
);

void build_weight_matrix(
    const uint8_t *A,
    double *S,
    int N, int Ne, int Ni,
    double max_exc_weight,
    double max_inh_weight
);

int make_output_filename(
    const char *adj_filename,
    double max_exc_weight,
    char *out_filename,
    int out_size, 
    int sim_number
);

int simulate_izhikevich(
    const uint8_t *A,
    int N,
    int SIM_TIME,
    double max_exc_weight,
    double max_inh_weight,
    double noise_max,
    const char *output_filename
);

int sticky_cylinder3D(double initial_x, double initial_y, double initial_z,
                      double *dx, double *dy, double *dz,
                      double L,
                      double *end_x, double *end_y, double *end_z,
                      double *X, double *Y, double *Z,
                      double *dendrites_diameters,
                      int i, int N_neurons,
                      uint8_t *AdjMatrix_flat,
                      int *trials, int *links,
                      FILE *axon_simulation,
                      int neuron_idx, int segment_idx);

double reflect_polar_angle(double angle_pol);
//************************************