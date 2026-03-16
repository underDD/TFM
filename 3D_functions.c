#include "3D_head.h"
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

void new_vector_segment_3D(double segment_length, double angle_azi, double angle_pol, double *dx, double *dy, double *dz){
    *dx = segment_length * cos(angle_azi) * sin(angle_pol);
    *dy = segment_length * sin(angle_azi) * sin(angle_pol);
    *dz = segment_length * cos(angle_pol);
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

bool new_axon_intersection_3D(double *x, double *y, double *z, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double *segment_vector_z, double L, double dx, double dy, double dz){

    double cx, cy, cz;

    cx = x[target_neuron] - *segment_vector_x;
    cy = y[target_neuron] - *segment_vector_y;
    cz = z[target_neuron] - *segment_vector_z;
    cx = min_image(cx, L);
    cy = min_image(cy, L);
    cz = min_image(cz, L);

    double f_raw = (cx * dx + cy * dy + cz * dz) / (dx * dx + dy * dy + dz * dz);
    double f = f_raw;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    // printf("f_raw: %f, f: %f\n", f_raw, f);

    double gx, gy, gz;

    gx = cx - f * dx;
    gy = cy - f * dy;
    gz = cz - f * dz;

    double distance = sqrt(gx*gx + gy*gy + gz*gz);

    if (distance <= dendrites_diameter[target_neuron]*0.5){
        return true;
    }

    return false;

}

bool new_axon_intersection_noPBC_3D(double cx, double cy, double cz, double dendrite_diameter,
                                double x0, double y0, double z0,
                                double dx, double dy, double dz){
    
    double R = dendrite_diameter / 2.0;
    
    // End of the segment
    double x1 = x0 + dx;
    double y1 = y0 + dy;
    double z1 = z0 + dz;

    // Vector of the segment
    double vx = x1 - x0; // = dx
    double vy = y1 - y0; // = dy
    double vz = z1 - z0; // = dz

    // Vector from the start of the segment to the center of the sphere
    double wx = cx - x0;
    double wy = cy - y0;
    double wz = cz - z0;

    // Length squared of the segment
    double vv = vx*vx + vy*vy + vz*vz;
    double t;

    // Segment is almost a point
    if (vv < 1e-15) {
        double ddx = cx - x0;
        double ddy = cy - y0;
        double ddz = cz - z0;
        return (ddx*ddx + ddy*ddy + ddz*ddz <= R*R);
    }

    // Projection of w onto v, normalized by the length of v
    t = (wx*vx + wy*vy + wz*vz) / vv;

    // t = 0 --> beginning of the segment, t = 1 --> end of the segment
    // t < 0 --> before the segment, t > 1 --> after the segment
    // 0 < t < 1 --> projection is on the segment
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    // Closest pont on the segment to the center of the sphere
    double px = x0 + t*vx;
    double py = y0 + t*vy;
    double pz = z0 + t*vz;

    // Distance from the closest point to the center of the sphere
    double ddx = cx - px;
    double ddy = cy - py;

    // Check if the distance is less than or equal to the radius of the sphere
    return (ddx*ddx + ddy*ddy + (cz - pz)*(cz - pz) <= R*R);

}

void sticky_walls3D(double initial_segment_vector_x, double initial_segment_vector_y, double initial_segment_vector_z, 
                    double *dx, double *dy, double *dz, double L,
                    double *end_segment_vector_x, double *end_segment_vector_y, double *end_segment_vector_z,
                    double *X, double *Y, double *Z, double *dendrites_diameters,
                    int i, int N_neurons, uint8_t *AdjMatrix_flat, double alpha,
                    int *trials, int *links, FILE *axon_simulation, int neuron_idx, int segment_idx){

    double xmin = -L/2.0, xmax = L/2.0;
    double ymin = -L/2.0, ymax = L/2.0;
    double zmin = -L/2.0, zmax = L/2.0;
    double eps = 1e-12;

    double x = initial_segment_vector_x;
    double y = initial_segment_vector_y;
    double z = initial_segment_vector_z;

    double total_len = sqrt((*dx)*(*dx) + (*dy)*(*dy) + (*dz)*(*dz));

    if (total_len < eps){
        *end_segment_vector_x = x;
        *end_segment_vector_y = y;
        *end_segment_vector_z = z;
        *dx = 0.0;
        *dy = 0.0;
        *dz = 0.0;
        return;
    }

    double ux = (*dx) / total_len;
    double uy = (*dy) / total_len;
    double uz = (*dz) / total_len;

    double remaining = total_len;

    int safety_counter = 0;
    int max_iter = 20;

    int subsegment_idx = 0;
    while (remaining > eps && safety_counter < max_iter){

        safety_counter++;

        double s_hit = DBL_MAX;
        int wall_type = 0; // 0 = no wall, 1 = x-wall, 2 = y-wall, 3 = z-wall

        // Distancia a caras x = const
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

        // Distancia a caras y = const
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

        // Distancia a caras z = const
        if (uz > eps){
            double s = (zmax - z) / uz;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 3;
            }
        }
        else if (uz < -eps){
            double s = (zmin - z) / uz;
            if (s >= -eps && s < s_hit){
                s_hit = s;
                wall_type = 3;
            }
        }

        double x_next, y_next, z_next;
        double step_len;

        // No choca antes de terminar la longitud restante
        if (s_hit == DBL_MAX || s_hit > remaining){
            x_next = x + remaining * ux;
            y_next = y + remaining * uy;
            z_next = z + remaining * uz;
            step_len = remaining;
            remaining = 0.0;
        }
        // Ya está prácticamente en una cara/arista/esquina
        else if (s_hit < eps){

            if (wall_type == 1) ux = 0.0;
            else if (wall_type == 2) uy = 0.0;
            else if (wall_type == 3) uz = 0.0;

            double norm_tan = sqrt(ux*ux + uy*uy + uz*uz);

            if (norm_tan < eps){
                break;
            }

            ux /= norm_tan;
            uy /= norm_tan;
            uz /= norm_tan;

            continue;
        }
        // Choca antes de terminar
        else{
            x_next = x + s_hit * ux;
            y_next = y + s_hit * uy;
            z_next = z + s_hit * uz;
            step_len = s_hit;
            remaining -= s_hit;
        }

        // Fijar exactamente en caras por estabilidad numérica
        if (fabs(x_next - xmax) < 1e-10) x_next = xmax;
        if (fabs(x_next - xmin) < 1e-10) x_next = xmin;
        if (fabs(y_next - ymax) < 1e-10) y_next = ymax;
        if (fabs(y_next - ymin) < 1e-10) y_next = ymin;
        if (fabs(z_next - zmax) < 1e-10) z_next = zmax;
        if (fabs(z_next - zmin) < 1e-10) z_next = zmin;

        // Subsegmento real recorrido
        double seg_dx = x_next - x;
        double seg_dy = y_next - y;
        double seg_dz = z_next - z;

        // Comprobar intersecciones sobre este subsegmento
        for (int j = 0; j < N_neurons; j++){
            if (j == i) continue;
            if (AdjMatrix_flat[i*N_neurons + j] != 0) continue;

            if (new_axon_intersection_noPBC_3D(
                    X[j], Y[j], Z[j], dendrites_diameters[j],
                    x, y, z,
                    seg_dx, seg_dy, seg_dz))
            {
                (*trials)++;
                if (randomInPR(0.0, 1.0) < alpha){
                    AdjMatrix_flat[i*N_neurons + j] = 1;
                    (*links)++;
                }
            }
        }

        if (axon_simulation != NULL){
            fprintf(axon_simulation, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
                    neuron_idx, segment_idx, subsegment_idx,
                    x, y, z,
                    x_next, y_next, z_next);
        }

        subsegment_idx++;

        // Avanzar al final del subsegmento
        x = x_next;
        y = y_next;
        z = z_next;

        // Si hemos tocado una cara y aún queda longitud, proyectar tangencialmente
        if (remaining > eps && s_hit != DBL_MAX){

            if (wall_type == 1) ux = 0.0;
            else if (wall_type == 2) uy = 0.0;
            else if (wall_type == 3) uz = 0.0;

            double norm_tan = sqrt(ux*ux + uy*uy + uz*uz);

            if (norm_tan < eps){
                break;
            }

            ux /= norm_tan;
            uy /= norm_tan;
            uz /= norm_tan;
        }
    }

    // Clamp final por seguridad
    if (x < xmin) x = xmin;
    if (x > xmax) x = xmax;
    if (y < ymin) y = ymin;
    if (y > ymax) y = ymax;
    if (z < zmin) z = zmin;
    if (z > zmax) z = zmax;

    *end_segment_vector_x = x;
    *end_segment_vector_y = y;
    *end_segment_vector_z = z;

    *dx = x - initial_segment_vector_x;
    *dy = y - initial_segment_vector_y;
    *dz = z - initial_segment_vector_z;
}