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

bool new_axon_intersection_3D(double *x, double *y, double *z, double *dendrites_diameter, int origin_neuron, int target_neuron, double *segment_vector_x, double *segment_vector_y, double *segment_vector_z, double L, double height, double dx, double dy, double dz){

    double cx, cy, cz;

    cx = x[target_neuron] - *segment_vector_x;
    cy = y[target_neuron] - *segment_vector_y;
    cz = z[target_neuron] - *segment_vector_z;
    cx = min_image(cx, L);
    cy = min_image(cy, L);
    cz = min_image(cz, height);

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
    double ddz = cz - pz;

    // Check if the distance is less than or equal to the radius of the sphere
    return (ddx*ddx + ddy*ddy + ddz*ddz <= R*R);

}

void sticky_walls3D(double initial_segment_vector_x, double initial_segment_vector_y, double initial_segment_vector_z, 
                    double *dx, double *dy, double *dz, double L, double height,
                    double *end_segment_vector_x, double *end_segment_vector_y, double *end_segment_vector_z,
                    double *X, double *Y, double *Z, double *dendrites_diameters,
                    int i, int N_neurons, uint8_t *AdjMatrix_flat,
                    int *trials, int *links, FILE *axon_simulation, int neuron_idx, int segment_idx){

    double xmin = -L/2.0, xmax = L/2.0;
    double ymin = -L/2.0, ymax = L/2.0;
    double zmin = -height/2.0, zmax = height/2.0;
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

            if (wall_type == 3){
                remaining = 0.0;
                break;
            }

            if (wall_type == 1) ux = 0.0;
            else if (wall_type == 2) uy = 0.0;

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
            // if (AdjMatrix_flat[i*N_neurons + j] != 0) continue;

            if (new_axon_intersection_noPBC_3D(
                    X[j], Y[j], Z[j], dendrites_diameters[j],
                    x, y, z,
                    seg_dx, seg_dy, seg_dz))
            {
                (*trials)++;
                
                    AdjMatrix_flat[i*N_neurons + j] += 1;
                    (*links)++;
                
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

        // Si toca la tapa superior o inferior, se corta el crecimiento
        if (remaining > eps && s_hit != DBL_MAX){

            if (wall_type == 3){
                remaining = 0.0;
                break;
            }

            if (wall_type == 1) ux = 0.0;
            else if (wall_type == 2) uy = 0.0;

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
    MKDIR("3D_dynamics");

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

    // Eliminar sufijo "_simX" si existe
    char *sim_ptr = strstr(stem, "_sim");
    if (sim_ptr) {
        *sim_ptr = '\0';
        }

    // Construir nombre final
    int written = snprintf(
        out_filename,
        out_size,
        "3D_dynamics/%s_exc%.3f_sim%d.txt",
        stem,
        max_exc_weight,
        sim_number
    );

    // printf("Output filename: %s\n", out_filename);

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

int sticky_cylinder3D(double initial_x, double initial_y, double initial_z,
                      double *dx, double *dy, double *dz,
                      double L, double height,
                      double *end_x, double *end_y, double *end_z,
                      double *X, double *Y, double *Z,
                      double *dendrites_diameters,
                      int i, int N_neurons,
                      uint8_t *AdjMatrix_flat,
                      int *trials, int *links,
                      FILE *axon_simulation,
                      int neuron_idx, int segment_idx)
{
    double eps = 1e-12;

    double R = L / 2.0;
    double zmin = -height / 2.0;
    double zmax =  height / 2.0;

    double x = initial_x;
    double y = initial_y;
    double z = initial_z;

    double total_len = sqrt((*dx)*(*dx) + (*dy)*(*dy) + (*dz)*(*dz));

    if (total_len < eps) {
        *end_x = x;
        *end_y = y;
        *end_z = z;
        *dx = *dy = *dz = 0.0;
        return 0;
    }

    double ux = (*dx) / total_len;
    double uy = (*dy) / total_len;
    double uz = (*dz) / total_len;

    double remaining = total_len;
    int subsegment_idx = 0;
    int stop_axon = 0;

    // Corrección si empieza ligeramente fuera del cilindro
    double r0 = sqrt(x*x + y*y);
    if (r0 > R) {
        x = x / r0 * R;
        y = y / r0 * R;
    }
    if (z < zmin) z = zmin;
    if (z > zmax) z = zmax;

    while (remaining > eps && !stop_axon) {

        double s_hit = remaining;
        int hit_type = 0; 
        // 0 = no hit, 1 = pared lateral, 2 = tapa/base

        // ======================================================
        // 1) Posible choque con tapa o base
        // ======================================================

        if (uz > eps) {
            double s = (zmax - z) / uz;
            if (s >= -eps && s < s_hit) {
                s_hit = s;
                hit_type = 2;
            }
        }
        else if (uz < -eps) {
            double s = (zmin - z) / uz;
            if (s >= -eps && s < s_hit) {
                s_hit = s;
                hit_type = 2;
            }
        }

        // ======================================================
        // 2) Posible choque con pared lateral: x² + y² = R²
        // ======================================================

        double a = ux*ux + uy*uy;
        double b = 2.0 * (x*ux + y*uy);
        double c = x*x + y*y - R*R;

        if (a > eps) {
            double disc = b*b - 4.0*a*c;

            if (disc >= 0.0) {
                double sqrt_disc = sqrt(disc);

                double s1 = (-b - sqrt_disc) / (2.0*a);
                double s2 = (-b + sqrt_disc) / (2.0*a);

                double s_wall = DBL_MAX;

                if (s1 > eps && s1 < s_wall) s_wall = s1;
                if (s2 > eps && s2 < s_wall) s_wall = s2;

                if (s_wall < s_hit && s_wall <= remaining + eps) {
                    s_hit = s_wall;
                    hit_type = 1;
                }
            }
        }

        // ======================================================
        // 3) Avanzar hasta el siguiente evento
        // ======================================================

        double x_next = x + s_hit * ux;
        double y_next = y + s_hit * uy;
        double z_next = z + s_hit * uz;

        // Correcciones numéricas
        double r_next = sqrt(x_next*x_next + y_next*y_next);
        if (r_next > R && r_next > eps) {
            x_next = x_next / r_next * R;
            y_next = y_next / r_next * R;
        }

        if (fabs(z_next - zmax) < 1e-10) z_next = zmax;
        if (fabs(z_next - zmin) < 1e-10) z_next = zmin;

        double seg_dx = x_next - x;
        double seg_dy = y_next - y;
        double seg_dz = z_next - z;

        // ======================================================
        // 4) Intersecciones en este subsegmento
        // ======================================================

        if (seg_dx*seg_dx + seg_dy*seg_dy + seg_dz*seg_dz > eps*eps) {

            for (int j = 0; j < N_neurons; j++) {

                if (j == i) continue;

                if (new_axon_intersection_noPBC_3D(
                        X[j], Y[j], Z[j],
                        dendrites_diameters[j],
                        x, y, z,
                        seg_dx, seg_dy, seg_dz
                    ))
                {
                    (*trials)++;
                    AdjMatrix_flat[i*N_neurons + j] += 1;
                    (*links)++;
                }
            }

            /*
            if (axon_simulation != NULL) {
                fprintf(axon_simulation,
                        "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
                        neuron_idx, segment_idx, subsegment_idx,
                        x, y, z,
                        x_next, y_next, z_next);
            }
            */

            subsegment_idx++;
        }

        remaining -= s_hit;

        x = x_next;
        y = y_next;
        z = z_next;

        // ======================================================
        // 5) Si toca tapa/base, se corta el axón
        // ======================================================

        if (hit_type == 2) {
            stop_axon = 1;
            break;
        }

        // ======================================================
        // 6) Si toca pared lateral, proyectar dirección sobre
        //    el plano tangente del cilindro
        // ======================================================

        if (remaining > eps && hit_type == 1) {

            double r = sqrt(x*x + y*y);

            if (r < eps) {
                break;
            }

            // Normal lateral del cilindro
            double nx = x / r;
            double ny = y / r;
            double nz = 0.0;

            // Quitar componente normal: u_tan = u - (u·n)n
            double dot = ux*nx + uy*ny + uz*nz;

            ux = ux - dot * nx;
            uy = uy - dot * ny;
            uz = uz - dot * nz;

            double norm_tan = sqrt(ux*ux + uy*uy + uz*uz);

            if (norm_tan < eps) {
                stop_axon = 1;
                break;
            }

            ux /= norm_tan;
            uy /= norm_tan;
            uz /= norm_tan;
        }
    }

    // Clamp final
    double rf = sqrt(x*x + y*y);
    if (rf > R && rf > eps) {
        x = x / rf * R;
        y = y / rf * R;
    }

    if (z < zmin) z = zmin;
    if (z > zmax) z = zmax;

    *end_x = x;
    *end_y = y;
    *end_z = z;

    *dx = x - initial_x;
    *dy = y - initial_y;
    *dz = z - initial_z;

    return stop_axon;
}

double reflect_polar_angle(double angle_pol) {
    angle_pol = fmod(angle_pol, 2.0 * PI);

    if (angle_pol < 0.0) {
        angle_pol += 2.0 * PI;
    }

    if (angle_pol > PI) {
        angle_pol = 2.0 * PI - angle_pol;
    }

    return angle_pol;
}