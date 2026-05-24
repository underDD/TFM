from pathlib import Path
from matplotlib.pylab import seed
import numpy as np
import matplotlib.pyplot as plt
from pypalettes import load_cmap 
import importlib
import functions_create_networks as fcn
importlib.reload(fcn)
from numba import njit
import sys

# fcn.compile_c_codes3D()

n_chunks = 5
chunk_id = int(sys.argv[1])

base_sigmas_all = np.linspace(0.04, 0.6, 50)
base_sigmas = np.array_split(base_sigmas_all, n_chunks)[chunk_id]

MAX_EXC_WEIGHT = np.linspace(2, 20.0, 50)
de = MAX_EXC_WEIGHT[1] - MAX_EXC_WEIGHT[0]

MAX_EXC_WEIGHT = np.arange(2, 30, de)

GNA_mean = []
GNA_var = []
GNA_fano = []


parameters3D = {
    "L": 3.0, # Length of the square domain (mm) or diameter of the circular domain (mm)
    "rho": 370, # Neuron density (neurons/mm^2)
    "soma_diameter": 0.015, # Diameter of the soma (mm)
    "d_mean": 0.15, # Dendritic tree diameter (mm) Gaussian distribution
    "d_sigma": 0.02, # Standard deviation of the dendritic tree diameter (mm)
    "l_mean": 0.30, # Mean axon length (mm) Rayleigh distribution
    "segment_length": 0.01, # Segment length for axon growth (mm)
    "sigma_pol": 0.1, # Standard deviation of the axon angle (radians)
    "sigma_azi": 0.1, # Standard deviation of the azimuthal
    "agg": 1, # Control parameter for the aggregation of neurons
    "n_centers": 400, # Number of centers for the aggregation of neurons
    "base_sigma": 0.8, # Base standard deviation for aggregation centers (mm)
    "BC_type": "SW",  # Boundary conditions type (PBC: Periodic Boundary Conditions, SW: Sticky Walls)
    "geometry": "cube" # Geometry of the domain (cube or cylinder)
}

alpha = 0.2

if parameters3D["geometry"] == "cylinder" and parameters3D["BC_type"] == "PBC":
    raise ValueError("Periodic Boundary Conditions (PBC) are not compatible with a cylindrical geometry. Please choose 'SW' for BC_type or change the geometry to 'cube'.")
if parameters3D["geometry"] == "cube":
    N = int(parameters3D["rho"] * parameters3D["L"]**3) # Total number of neurons
elif parameters3D["geometry"] == "cylinder":
    R = parameters3D["L"] / np.sqrt(np.pi)
    print(f"Effective radius of the cylinder domain: {R:.2f} mm")
    N = int(parameters3D["rho"] * np.pi * R**2 * parameters3D["L"]) # Total number of neurons, same as in the cube case for the same L, but we keep it here for clarity

print("Geometry:", parameters3D["geometry"])
print(f"Domain size (L): {parameters3D['L']} mm")
print(f"Neuron density (rho): {parameters3D['rho']} neurons/mm^3")
print(f"Number of neurons: {N}")

parametersDynamics = {
    "max_exc_weight": 8,
    "max_inh_weight": 0.5,
    "noise_max": 3,
    "sim_time": 10000,
    "seed": 42
}

new_sim = True
sims = 1

# output_file = f"Gini_things/3D_GNAmap_{parameters3D['geometry']}_{parameters3D['BC_type']}_l{parameters3D['l_mean']:.2f}_sim{sim}.txt"
# open(output_file, "w").close()  # Limpiar el archivo antes de escribir
# with open(output_file, "w") as f:
#     f.write("base_sigma\tMAX_EXC_WEIGHT\tGini\tmean_peaks_height\tstd_peaks_height\tfano_factor\n")

Gini = 0
cell_size = 0.25 # (mm)
n_grid = int(np.round(parameters3D["L"] / cell_size))


for sim in range(sims):
    output_file = f"Gini_things/3D_GNAmap_{parameters3D['geometry']}_{parameters3D['BC_type']}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_sim{sim}_full_chunk{chunk_id}.txt"
    open(output_file, "w").close()  # Limpiar el archivo antes de escribir
    with open(output_file, "w") as f:
        f.write("base_sigma\tMAX_EXC_WEIGHT\tGini\tmean_peaks_height\tstd_peaks_height\tfano_factor\n")
    for s in base_sigmas:

        parameters3D["base_sigma"] = s
        filename_neuron_params = fcn.run_3D_tune_positions(parameters3D, new_sim=new_sim)
        filename_network_tries = fcn.create_network3D(parameters3D, new_sim=new_sim)
        # X, Y, Soma_Diameter, Dendrite_Diameter, Axon_Length =fcn.load_neuron_params(filename_neuron_params)
        # Gini = fcn.compute_spatial_gini(X, Y, parameters3D["L"], geometry=parameters3D["geometry"], n_grid=n_grid)[0]

        A, filename_adj = fcn.modify_A_alpha3D(
            fcn.load_A(filename_network_tries),
            alpha,
            parameters3D,
            new_sim=new_sim,
            sim = sim
        )

        N = A.shape[0]

        for e in MAX_EXC_WEIGHT:

            parametersDynamics["max_exc_weight"] = e
            parametersDynamics["out_filename"] = filename_adj

            filename_dynamics = fcn.run_dynamics_3D(
                parametersDynamics,
                parameters3D,
                alpha,
                sim_number=sim,
                new_sim=new_sim, 
                show_prints_control=False
            )

            firings_t, firings_i = fcn.load_firings(filename_dynamics)

            assert firings_i.max() <= N

            GNA = fcn.compute_GNA_from_raster(
                    firings_t=firings_t,
                    firings_i=firings_i,
                    N=N,
                    SIM_TIME=parametersDynamics["sim_time"],
                    W=25
                )
            peaks, peak_heights, mean_peaks_height, std_peaks_height = fcn.find_GNA_peaks(GNA) 
            fano_factor = std_peaks_height**2 / mean_peaks_height if mean_peaks_height > 0 else 0.0

            GNA_mean.append(GNA.mean())
            GNA_var.append(GNA.var())
            GNA_fano.append(fano_factor)

            with open(output_file, "a") as f:
                f.write(f"{s:.4f}\t{e:.2f}\t{Gini:.6f}\t{mean_peaks_height:.6f}\t{std_peaks_height:.6f}\t{fano_factor:.6f}\n")
            
            print(f"base_sigma: {s:.4f}, MAX_EXC_WEIGHT: {e:.2f}, Gini: {Gini:.6f}, mean_peaks_height: {mean_peaks_height:.6f}, std_peaks_height: {std_peaks_height:.6f}, fano_factor: {fano_factor:.6f}")
    print(f"Simulation {sim+1}/{sims} completed, saved in {output_file}", end="\n")

# Le cuesta 10h
# Mirar distribuciones de kIn con las dos condiciones periodicas de contorno
# Plotear las neuronas y cambiar e tamaño del nodo en funcion del grado para cada condicones, hacer smoothing de los colores
# Poner en el mapa ruster en extreos e intemedios
# Grupo de neuronas inactivas
# Hacer el cultivo cilindrico con SW en las paredes y en las bases que se corte el crecimiento. 
# Mirar el small_world, richclub