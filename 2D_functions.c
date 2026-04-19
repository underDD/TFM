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
    @param trials: pointer to count the number of intersection trials
    @param links: pointer to count the number of successful connections formed
*/
void sticky_walls2D(double initial_segment_vector_x, double initial_segment_vector_y, double *dx, double *dy,
                       double L, double *end_segment_vector_x, double *end_segment_vector_y, double *X, double *Y,
                       double *dendrites_diameters, int i, int N_neurons, uint8_t *AdjMatrix_flat,
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
            
                    AdjMatrix_flat[i*N_neurons + j] = 1;
                    (*links)++;
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

int get_matrix_size(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening adjacency file");
        return -1;
    }

    int rows = 0;
    char line[MAX_STRING_LENGTH];

    while (fgets(line, sizeof(line), file)) {
        if (strlen(line) > 1) {
            rows++;
        }
    }

    fclose(file);
    return rows;
}

int load_adjacency_matrix(const char *filename, uint8_t *AdjMatrix_flat, int N) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening adjacency file");
        return -1;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int val;
            if (fscanf(file, "%d", &val) != 1) {
                fprintf(stderr, "Error reading adjacency matrix at row %d col %d\n", i, j);
                fclose(file);
                return -1;
            }
            AdjMatrix_flat[i*N + j] = (uint8_t)val;
        }
    }

    fclose(file);
    return 0;
}

/*
    Initialize the Izhikevich parameters for excitatory and inhibitory neurons with some variability.
    Ne: number of excitatory neurons
    Ni: number of inhibitory neurons
    a, b, c, d: arrays to store the parameters for each neuron
*/
void init_izhikevich_parameters(
    int N, int Ne, int Ni,
    double *a, double *b, double *c, double *d
) {
    for (int i = 0; i < Ne; i++) { // Para las excitatorias, a es 0.02, b es 0.2, c varía entre -65 y -50, d varía entre 2 y 8
        double re = randomInPR(0.0, 1.0);
        a[i] = 0.02;
        b[i] = 0.2;
        c[i] = -65.0 + 15.0 * re * re;
        d[i] = 8.0 - 6.0 * re * re;
    }

    for (int i = 0; i < Ni; i++) { // Para las inhibitorias, a varía entre 0.02 y 0.1, b varía entre 0.25 y 0.2, c es -65, d es 2
        double ri = randomInPR(0.0, 1.0);
        int idx = Ne + i;
        a[idx] = 0.02 + 0.08 * ri;
        b[idx] = 0.25 - 0.05 * ri;
        c[idx] = -65.0;
        d[idx] = 2.0;
    }
}

/*  
    Build the synaptic weight matrix S based on the adjacency matrix A and the maximum excitatory and inhibitory weights.
    A: adjacency matrix (flattened)
    S: weight matrix to be filled (flattened)
    N: total number of neurons
    Ne: number of excitatory neurons
    Ni: number of inhibitory neurons
    max_exc_weight: maximum weight for excitatory connections
    max_inh_weight: maximum weight for inhibitory connections
*/
void build_weight_matrix(
    const uint8_t *A,
    double *S,
    int N, int Ne, int Ni,
    double max_exc_weight,
    double max_inh_weight
) {
    for (int i = 0; i < N; i++) { // Para cada neurona i, asignamos pesos a sus conexiones salientes
        for (int j = 0; j < N; j++) {
            double w;
            
            // Las primeras Ne neuronas son excitatorias, las últimas Ni son inhibitorias
            if (j < Ne) { // Si la neurona j es excitatoria, el peso es positivo y se asigna un valor aleatorio entre 0 y max_exc_weight
                w = max_exc_weight * randomInPR(0.0, 1.0);
            } else { // Si la neurona j es inhibitoria, el peso es negativo y se asigna un valor aleatorio entre -max_inh_weight y 0
                w = -max_inh_weight * randomInPR(0.0, 1.0);
            }

            S[i*N + j] = ((double)A[i*N + j]) * w;
        }
    }
}

#include <stdio.h>
#include <string.h>

int make_output_filename(
    const char *adj_filename,
    double max_exc_weight,
    char *out_filename,
    int out_size, 
    int sim_number
) {
    MKDIR("2D_dynamics");

    // Obtener nombre base (sin ruta)
    const char *base = strrchr(adj_filename, '/');
#ifdef _WIN32
    const char *base2 = strrchr(adj_filename, '\\');
    if (base2 && (!base || base2 > base)) base = base2;
#endif
    base = base ? base + 1 : adj_filename;

    // Copiar a buffer local
    char stem[512];
    strncpy(stem, base, sizeof(stem) - 1);
    stem[sizeof(stem) - 1] = '\0';

    // Quitar extensión (.txt)
    char *dot = strrchr(stem, '.');
    if (dot) *dot = '\0';

    // Sustituir "adjacency_matrix" por "dynamics"
    char *p = strstr(stem, "adjacency_matrix");
    if (p) {
        memmove(p + strlen("dynamics"),
                p + strlen("adjacency_matrix"),
                strlen(p + strlen("adjacency_matrix")) + 1);
        memcpy(p, "dynamics", strlen("dynamics"));
    }

    // Construir nombre final
    int written = snprintf(
        out_filename,
        out_size,
        "2D_dynamics/%s_exc%.3f_sim%d.txt",
        stem,
        max_exc_weight,
        sim_number
    );

    if (written < 0 || written >= out_size) {
        fprintf(stderr, "Output filename too long.\n");
        return -1;
    }

    return 0;
}

/*
    Simulate the Izhikevich neuron model dynamics based on the adjacency matrix A and save the spike times to an output file.
    A: adjacency matrix (flattened)
    N: total number of neurons
    SIM_TIME: total simulation time in ms
    max_exc_weight: maximum excitatory weight (used for building the weight matrix)
    max_inh_weight: maximum inhibitory weight (used for building the weight matrix)
    noise_max: maximum noise amplitude for excitatory neurons
    output_filename: name of the file to save spike times
*/
int simulate_izhikevich(
    const uint8_t *A,
    int N,
    int SIM_TIME,
    double max_exc_weight,
    double max_inh_weight,
    double noise_max,
    const char *output_filename
) {

    // 80% excitatory, 20% inhibitory
    int Ne = (int)(0.8 * N);
    int Ni = N - Ne;

    // Allocate memory for neuron parameters and state variables
    double *a = (double *)malloc(N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *c = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));
    double *v = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));
    double *I = (double *)malloc(N * sizeof(double));
    double *S = (double *)malloc(N * N * sizeof(double));
    int *fired = (int *)malloc(N * sizeof(int));

    // Esto es para verificar que la memoria se ha asignado correctamente
    if (!a || !b || !c || !d || !v || !u || !I || !S || !fired) {
        fprintf(stderr, "Memory allocation error in simulate_izhikevich.\n");
        free(a); free(b); free(c); free(d); free(v); free(u); free(I); free(S); free(fired);
        return -1;
    }

    // Inicializa los parámetros de las neuronas a, b, c, d y construye la matriz sináptica S a partir de la matriz de adyacencia A
    init_izhikevich_parameters(N, Ne, Ni, a, b, c, d);
    build_weight_matrix(A, S, N, Ne, Ni, max_exc_weight, max_inh_weight);

    // Inicializa las variables de estado v y u
    for (int i = 0; i < N; i++) {
        v[i] = -65.0;
        u[i] = b[i] * v[i];
    }

    FILE *out = fopen(output_filename, "w");
    if (!out) {
        perror("Error opening output file");
        free(a); free(b); free(c); free(d); free(v); free(u); free(I); free(S); free(fired);
        return -1;
    }

    fprintf(out, "# time\tneuron\n");

    // long long para contar el total de spikes, ya que puede ser un número muy grande
    long long total_spikes = 0;

    // Simulación del modelo de Izhikevich
    for (int t = 1; t <= SIM_TIME; t++) {

        // Ruido gaussiano para las neuronas excitatorias, ruido más pequeño para las inhibitorias
        for (int i = 0; i < Ne; i++) {
            I[i] = noise_max * box_muller();
        }
        for (int i = Ne; i < N; i++) {
            I[i] = 2.0 * box_muller();
        }

        // Detectar qué neuronas han disparado (v >= 30) y almacenar sus índices en el array fired
        int n_fired = 0;
        for (int i = 0; i < N; i++) {
            if (v[i] >= 30.0) {
                fired[n_fired++] = i;
            }
        }

        // Para cada neurona que ha disparado se guarda en el fichero de salida y se restea su valor de v y u
        if (n_fired > 0) {
            for (int k = 0; k < n_fired; k++) {
                int idx = fired[k];

                fprintf(out, "%d\t%d\n", t, idx + 1); // Guardamos el tiempo y el índice de la neurona (sumamos 1 para que empiece en 1 en lugar de 0)
                total_spikes++;

                v[idx] = c[idx]; // Reset del potencial de membrana
                u[idx] += d[idx]; // Reset de la variable de recuperación
            }

            // Para cada neurona i se suman todas las señales que llegan de las neuronas que han disparado.
            for (int i = 0; i < N; i++) {
                double sum_input = 0.0; // Para acumular lo que llega a i
                for (int k = 0; k < n_fired; k++) { // Bucle sobre las neuronas que han disparado
                    int j = fired[k]; // Índice de la neurona que ha disparado
                    sum_input += S[i*N + j]; // S[i*N + j] es el peso de la conexión de j a i, si no hay conexión es 0, si j es excitatoria es positivo, si j es inhibitoria es negativo
                }
                I[i] += sum_input; // Se añade la suma de las entradas a la corriente I[i] que ya tiene el ruido
            }
        }

        // Integrar las ecuaciones
        for (int i = 0; i < N; i++) {
            v[i] += 0.5 * (0.04*v[i]*v[i] + 5.0*v[i] + 140.0 - u[i] + I[i]);
            v[i] += 0.5 * (0.04*v[i]*v[i] + 5.0*v[i] + 140.0 - u[i] + I[i]);
            u[i] += a[i] * (b[i]*v[i] - u[i]);
        }
    }

    fclose(out);

    printf("Number of neurons: %d\n", N);
    printf("Excitatory neurons: %d\n", Ne);
    printf("Inhibitory neurons: %d\n", Ni);
    printf("Total spikes: %lld\n", total_spikes);
    printf("Dynamics saved in: %s\n", output_filename);

    free(a); free(b); free(c); free(d); free(v); free(u); free(I); free(S); free(fired);
    return 0;
}