#include "head.h"
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

/*
    Generate a random number between min and max
*/
float randomIn(double min, double max) {
    return min + (max - min) * fran;
}

/*
    Initialize the random number generator
    SEMILLA: seed
*/
void  ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }

    ind_ran=ig1=ig2=ig3=0;
}

/*
    Generate a random number between 0 and 1 using the Parisi-Rapuano method
*/
float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;

    return r;
}

/*
    Generate a random number between min and max
*/
float randomInPR(double min, double max) {
    return min + (max - min) * Random();
}


/*
    Generate a random number using the Box-Muller method
*/
double box_muller(){

    double r1, r2, s, u1;
    
    // Generate two uniform random numbers between 0 and 1
    r1 = Random();
    r2 = Random();

    // Apply the Box-Muller transform
    s = sqrt(-2.0 * log(r1));
    u1 = s * cos(2.0 * PI * r2);
    //*u2 = s * sin(2.0 * M_PI * r2);
    return u1;

}

double sqrt_distance(double x, double y){
    return sqrt(x*x + y*y);
}

double scalar_product(double x1, double y1, double x2, double y2){
    return x1*x2 + y1*y2;
}

void PBC(double *x, double L){
    if (L <= 0.0) return;
    double r = fmod(*x + L/2.0, L);   // desplazar para obtener resto en torno a 0..L
    if (r < 0.0) r += L;              // fmod puede devolver negativo
    *x = r - L/2.0;                   // volver al intervalo [-L/2, L/2)
}

double inverse_cumulative_rayleigh(double p, double sigma){
    // printf("p: %f, sigma: %f\n", p, sigma);
    return sigma * sqrt(-2.0 * log(1.0 - p));
}

void new_vector_segment(double segment_length, double angle, double *dx, double *dy){
    *dx = segment_length * cos(angle);
    *dy = segment_length * sin(angle);
}

double calculate_f(double *x, double *y, int neuron, double segment_vector_x, double segment_vector_y){

    double cs = scalar_product(x[neuron], y[neuron], segment_vector_x, segment_vector_y);
    double ss = scalar_product(segment_vector_x, segment_vector_y, segment_vector_x, segment_vector_y);
    
    return cs / ss;
}

bool new_axon_intersection(double *x, double *y, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double segment_length, double sigma_axon_angle, double L){

    double dx,dy;
    double angle;
    double end_segment_x, end_segment_y;

    angle = box_muller() * sigma_axon_angle; // Standard deviation of sigma_axon_angle radians
    new_vector_segment(segment_length, angle, &dx, &dy);

    double f = calculate_f(x, y, target_neuron, dx, dy);
    double gx, gy;
    
    end_segment_x = *segment_vector_x + dx;
    end_segment_y = *segment_vector_y + dy;
    PBC(&end_segment_x, L);
    PBC(&end_segment_y, L);

    gx = *segment_vector_x + f * dx;
    gy = *segment_vector_y + f * dy;

    double distance = sqrt_distance(gx-x[target_neuron], gy-y[target_neuron]);

    if (distance <= dendrites_diameter[target_neuron]/2.0)
        return true;
    
    return false;
}