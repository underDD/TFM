#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <stdbool.h>
#include <float.h>


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

//**********functions.c*************  

float randomIn(double min, double max);
float Random(void);
void ini_ran(int SEMILLA);
float randomInPR(double min, double max);
double box_muller(void);
double sqrt_distance(double x, double y);
double scalar_product(double x1, double y1, double x2, double y2);
void PBC(double *x, double L);
double inverse_cumulative_rayleigh(double p, double sigma);
void new_vector_segment(double segment_length, double angle, double *dx, double *dy);
double calculate_f(double *x, double *y, int neuron, double segment_vector_x, double segment_vector_y);
bool new_axon_intersection(double *x, double *y, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double segment_length, double sigma_axon_angle, double L);

//************************************