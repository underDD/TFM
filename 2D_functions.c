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

/*
    Check if a line segment (x0,y0) -> (x0+dx, y0+dy) intersects with a circle centered at (cx,cy) with radius = dendrite_diameter/2
    @param cx, cy: center of the circle (dendrite)
    @param dendrite_diameter: diameter of the dendrite (circle)
    @param x0, y0: starting point of the line segment (axon)
    @param dx, dy: vector of the line segment (axon)
    @return: true if intersects, false otherwise
*/
bool new_axon_intersection_noPBC2D(double cx, double cy, double dendrite_diameter,
                                double x0, double y0,
                                double dx, double dy)
{
    double R = dendrite_diameter / 2.0;

    // End of the segment
    double x1 = x0 + dx;
    double y1 = y0 + dy;

    // Vector of the segment
    double vx = x1 - x0; // = dx
    double vy = y1 - y0; // = dy

    // Vector from the start of the segment to the center of the circle
    double wx = cx - x0;
    double wy = cy - y0;

    // Length squared of the segment vector
    double vv = vx*vx + vy*vy;
    double t;

    // Segment is almost a point
    if (vv < 1e-15){
        double ddx = cx - x0;
        double ddy = cy - y0;
        return (ddx*ddx + ddy*ddy <= R*R);
    }

    // Projection of w onto v, normalized by the length of v
    t = (wx*vx + wy*vy) / vv;

    // t = 0 --> beginning of the segment, t = 1 --> end of the segment
    // t < 0 --> before the segment, t > 1 --> after the segment
    // 0 < t < 1 --> projection is on the segment
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    // Closest point on the segment to the center of the circle
    double px = x0 + t*vx;
    double py = y0 + t*vy;

    // Distance from the closest point to the center of the circle
    double ddx = cx - px;
    double ddy = cy - py;

    // Check if the distance is less than or equal to the radius of the circle
    return (ddx*ddx + ddy*ddy <= R*R);
}

/*
    Simulate the growth of an axon segment with "sticky walls" boundary conditions.
    The axon grows in a straight line defined by (dx, dy) until it hits a wall. Upon hitting a wall, it "sticks" and continues growing along the wall until it can grow again without immediately hitting the wall.
    During the growth, it checks for intersections with dendrites of other neurons and updates the adjacency matrix accordingly.
    @param initial_segment_vector_x, initial_segment_vector_y: starting point of the axon segment
    @param dx, dy: initial growth vector of the axon segment (will be modified to reflect actual growth)
    @param L: size of the box (for boundary conditions)
    @param end_segment_vector_x, end_segment_vector_y: output parameters to store the final position of the axon segment after growth
    @param X, Y: arrays of neuron positions
    @param dendrites_diameters: array of dendrite diameters for each neuron
    @param i: index of the current neuron whose axon is growing
    @param N_neurons: total number of neurons
    @param AdjMatrix_flat: flattened adjacency matrix to update with new connections
    @param alpha: probability of forming a connection upon intersection
    @param trials: pointer to count the number of intersection trials
    @param links: pointer to count the number of successful connections formed
*/
void sticky_walls2D(double initial_segment_vector_x, double initial_segment_vector_y, double *dx, double *dy,
                       double L, double *end_segment_vector_x, double *end_segment_vector_y, double *X, double *Y,
                       double *dendrites_diameters, int i, int N_neurons, uint8_t *AdjMatrix_flat, double alpha,
                       int *trials, int *links, FILE *axon_simulation, int neuron_idx, int segment_idx){

    double xmin = -L/2.0, xmax = L/2.0;
    double ymin = -L/2.0, ymax = L/2.0;
    double eps = 1e-12;

    double x = initial_segment_vector_x;
    double y = initial_segment_vector_y;
    double total_len = sqrt(*dx*(*dx) + *dy*(*dy));
    
    // Normalized lenghts
    double ux = *dx / total_len; // The new position is x_new = x_old + ux * step_len, where step_len is how much we can grow until hitting a wall or exhausting the remaining length
    double uy = *dy / total_len; // y_new = y_old + uy * step_len
    double remaining = total_len; // remaining length to grow

    int subsegment_idx = 0;
    while (remaining > eps){

        double s_hit = DBL_MAX; // Start like the wall is in infinite
        int wall_type = 0; // 0 = no wall, 1 = vertical wall, 2 = horizontal wall

        // This four conditions are to check which is the closest wall

        if (ux > eps){ // Moving to the right, right wall case
            double s = (xmax - x)/ux; // how much we can grow until hitting the right wall
            if (s < s_hit){ // we hit the wall before exhausting the remaining length
                s_hit = s;
                wall_type = 1;
            }
        }
        else if (ux < -eps){ // Moving to the left, left wall case
            double s = (xmin - x)/ux; // how much we can grow until hitting the left wall
            if (s < s_hit){ // we hit the wall before exhausting the remaining length
                s_hit = s;
                wall_type = 1;
            }
        }

        if (uy > eps){ // Moving upwards, top wall case
            double s = (ymax - y)/uy; // how much we can grow until hitting the top wall
            if (s < s_hit){ // we hit the wall before exhausting the remaining length
                s_hit = s;
                wall_type = 2;
            }
        }
        else if (uy < -eps){ // Moving downwards, bottom wall case
            double s = (ymin - y)/uy; // how much we can grow until hitting the bottom wall
            if (s < s_hit){ // we hit the wall before exhausting the remaining length
                s_hit = s;
                wall_type = 2;
            }
        }

        double x_next, y_next;
        double step_len;

        if(s_hit > remaining){ // Case that the wall is far
            step_len = remaining;
        }
        else if(s_hit < eps){ // It is in the wall practicamente 
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
        else{ // It hits the wall in this step
            step_len = s_hit;
        }

        
        // The next step is just the sum of the step_lenght in the (ux, uy) direction
        x_next = x + step_len*ux;
        y_next = y + step_len*uy;
        remaining -= step_len;
        
        if (fabs(x_next -xmax) < eps) x_next = xmax;
        if (fabs(x_next -xmin) < eps) x_next = xmin;
        if (fabs(y_next -ymax) < eps) y_next = ymax;
        if (fabs(y_next -ymin) < eps) y_next = ymin;

        // Now check intersections in this new segment

        double seg_dx = x_next - x; // x = initial_x
        double seg_dy = y_next - y; // y = initial_y

        // This is to check for intersections to this new segment
        for (int j = 0; j < N_neurons; j++){
            // Check intersection with dendrite j
            if (j==i) continue;
            if (AdjMatrix_flat[i*N_neurons + j]!= 0) continue;

            if(new_axon_intersection_noPBC2D(X[i], Y[i], dendrites_diameters[j], x, y, seg_dx, seg_dy)){
                (*trials)++;
                if (randomInPR(0,1)<alpha){
                    AdjMatrix_flat[i*N_neurons + j] = 1;
                    (*links)++;
                }
            }
        }

        if (axon_simulation != NULL){
            fprintf(axon_simulation, "%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
                    neuron_idx, segment_idx, subsegment_idx,
                    x, y,
                    x_next, y_next);
        }

        subsegment_idx++;

        // Go to the end of this segment
        x = x_next; 
        y = y_next; 

        // If we just hit the wall, we need to reorient the growth vector to be parallel to the wall for the remaining length
        if (remaining > eps && s_hit != DBL_MAX && step_len == s_hit){
            if (wall_type == 1){
                ux = 0.0;
                if (uy > eps) uy = 1.0;
                else if (uy < -eps) uy = -1.0;
            }
            else if (wall_type == 2){
                uy = 0.0;
                if (ux > eps) ux = 1.0;
                else if (ux < -eps) ux = -1.0;
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