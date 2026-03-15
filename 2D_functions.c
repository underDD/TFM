#include "2D_head.h"
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


double inverse_cumulative_rayleigh(double p, double sigma){
    // printf("p: %f, sigma: %f\n", p, sigma);
    return sigma * sqrt(-2.0 * log(1.0 - p));
}

void new_vector_segment(double segment_length, double angle, double *dx, double *dy){
    *dx = segment_length * cos(angle);
    *dy = segment_length * sin(angle);
}

void PBC(double *x, double L){
    double r = fmod(*x + L/2.0, L);   // desplazar para obtener resto en torno a 0..L
    if (r < 0.0) r += L;              // fmod puede devolver negativo
    *x = r - L/2.0;                   // volver al intervalo [-L/2, L/2)
}

/*
    min_image: Calculate the minimum image distance considering periodic boundary conditions
    d: distance component (x or y)
    L: box length
*/
double min_image(double d, double L){
    if (d >  0.5*L) d -= L;
    if (d < -0.5*L) d += L;
    return d;
}

bool new_axon_intersection(double *x, double *y, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double L, double dx, double dy){

    double cx, cy;

    cx = x[target_neuron] - *segment_vector_x;
    cy = y[target_neuron] - *segment_vector_y;
    cx = min_image(cx, L);
    cy = min_image(cy, L);

    double f_raw = (cx * dx + cy * dy) / (dx * dx + dy * dy);
    double f = f_raw;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    // printf("f_raw: %f, f: %f\n", f_raw, f);

    double gx, gy;

    gx = cx - f * dx;
    gy = cy - f * dy;

    double distance = sqrt_distance(gx, gy);

    if (distance <= dendrites_diameter[target_neuron]*0.5){
        return true;
    }

    return false;

}

#include <math.h>
#include <float.h>
#include <stdint.h>

int new_axon_intersection_noPBC(double cx, double cy, double dendrite_diameter,
                                double x0, double y0,
                                double dx, double dy)
{
    double R = dendrite_diameter / 2.0;

    double x1 = x0 + dx;
    double y1 = y0 + dy;

    double vx = x1 - x0;
    double vy = y1 - y0;

    double wx = cx - x0;
    double wy = cy - y0;

    double vv = vx*vx + vy*vy;
    double t;

    if (vv < 1e-15){
        double ddx = cx - x0;
        double ddy = cy - y0;
        return (ddx*ddx + ddy*ddy <= R*R);
    }

    t = (wx*vx + wy*vy) / vv;

    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    double px = x0 + t*vx;
    double py = y0 + t*vy;

    double ddx = cx - px;
    double ddy = cy - py;

    return (ddx*ddx + ddy*ddy <= R*R);
}

void sticky_walls(double initial_segment_vector_x, double initial_segment_vector_y,
                  double *dx, double *dy, double L,
                  double *end_segment_vector_x, double *end_segment_vector_y,
                  double *X, double *Y, double *dendrites_diameters,
                  int i, int N_neurons,
                  uint8_t *AdjMatrix_flat, double alpha,
                  int *trials, int *links)
{
    double xmin = -L/2.0, xmax = L/2.0;
    double ymin = -L/2.0, ymax = L/2.0;
    double eps = 1e-12;

    double x = initial_segment_vector_x;
    double y = initial_segment_vector_y;

    double total_len = sqrt((*dx)*(*dx) + (*dy)*(*dy));

    if (total_len < eps){
        *end_segment_vector_x = x;
        *end_segment_vector_y = y;
        *dx = 0.0;
        *dy = 0.0;
        return;
    }

    double ux = (*dx) / total_len;
    double uy = (*dy) / total_len;
    double remaining = total_len;

    int safety_counter = 0;
    int max_iter = 20;

    while (remaining > eps && safety_counter < max_iter){

        safety_counter++;

        double s_hit = DBL_MAX;
        int wall_type = 0;   // 1 = vertical, 2 = horizontal

        // Distancia hasta pared vertical
        if (ux > eps){
            double s = (xmax - x) / ux;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 1;
            }
        }
        else if (ux < -eps){
            double s = (xmin - x) / ux;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 1;
            }
        }

        // Distancia hasta pared horizontal
        if (uy > eps){
            double s = (ymax - y) / uy;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 2;
            }
        }
        else if (uy < -eps){
            double s = (ymin - y) / uy;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 2;
            }
        }

        double x_next, y_next;
        double step_len;

        // Caso 1: no choca antes de agotar la longitud restante
        if (s_hit == DBL_MAX || s_hit > remaining){
            x_next = x + remaining * ux;
            y_next = y + remaining * uy;
            step_len = remaining;
            remaining = 0.0;
        }
        // Caso 2: está prácticamente ya en pared
        else if (s_hit < eps){
            if (wall_type == 1){
                ux = 0.0;
                if (uy > eps) uy = 1.0;
                else if (uy < -eps) uy = -1.0;
                else break;
            }
            else if (wall_type == 2){
                uy = 0.0;
                if (ux > eps) ux = 1.0;
                else if (ux < -eps) ux = -1.0;
                else break;
            }
            continue;
        }
        // Caso 3: choca con pared antes de terminar este trozo
        else{
            x_next = x + s_hit * ux;
            y_next = y + s_hit * uy;
            step_len = s_hit;
            remaining -= s_hit;
        }

        // Clamp fino por estabilidad
        if (fabs(x_next - xmax) < 1e-10) x_next = xmax;
        if (fabs(x_next - xmin) < 1e-10) x_next = xmin;
        if (fabs(y_next - ymax) < 1e-10) y_next = ymax;
        if (fabs(y_next - ymin) < 1e-10) y_next = ymin;

        // =========================================================
        // COMPROBAR INTERSECCIONES EN ESTE SUBSEGMENTO REAL: (x,y) -> (x_next,y_next)
        // =========================================================
        double seg_dx = x_next - x;
        double seg_dy = y_next - y;

        for (int j = 0; j < N_neurons; j++){
            if (j == i) continue;
            if (AdjMatrix_flat[i*N_neurons + j] != 0) continue;

            if (new_axon_intersection_noPBC(
                    X[j], Y[j], dendrites_diameters[j],
                    x, y,
                    seg_dx, seg_dy))
            {
                (*trials)++;
                if (randomInPR(0.0, 1.0) < alpha){
                    AdjMatrix_flat[i*N_neurons + j] = 1;
                    (*links)++;
                }
            }
        }

        // Avanzar realmente al final de este subsegmento
        x = x_next;
        y = y_next;

        // Si justo acabamos de tocar pared, reorientar el resto pegado a la pared
        if (remaining > eps && s_hit != DBL_MAX && step_len == s_hit){
            if (wall_type == 1){
                ux = 0.0;
                if (uy > eps) uy = 1.0;
                else if (uy < -eps) uy = -1.0;
                else break;
            }
            else if (wall_type == 2){
                uy = 0.0;
                if (ux > eps) ux = 1.0;
                else if (ux < -eps) ux = -1.0;
                else break;
            }
        }
    }

    if (x < xmin) x = xmin;
    if (x > xmax) x = xmax;
    if (y < ymin) y = ymin;
    if (y > ymax) y = ymax;

    *end_segment_vector_x = x;
    *end_segment_vector_y = y;

    *dx = x - initial_segment_vector_x;
    *dy = y - initial_segment_vector_y;
}