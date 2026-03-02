#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>


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
//************************************