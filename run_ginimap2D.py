from pathlib import Path
import sys
from matplotlib.pylab import seed
import numpy as np
import matplotlib.pyplot as plt
from pypalettes import load_cmap 
import importlib
import functions_create_networks as fcn
importlib.reload(fcn)
from numba import njit

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# fcn.compile_c_codes()

calibration_filename = "figures_TFM/data_for_figures/gini2D_n_grid13_calibration.txt"
data = np.loadtxt(calibration_filename, skiprows=1)
base_sigmas_calibration = data[:,1]
gini_calibration = data[:, 0]

n_chunks = 1
chunk_id = int(sys.argv[1])

# base_sigmas = np.array_split(base_sigmas_calibration, n_chunks)[chunk_id]
# print(f"Base sigmas to process in this chunk: {base_sigmas}")

base_sigmas = np.linspace(0.03, 0.7, 50)
base_sigmas = np.array_split(base_sigmas, n_chunks)[chunk_id]
MAX_EXC_WEIGHT = np.linspace(5, 30.0, 50)
# de = MAX_EXC_WEIGHT[1] - MAX_EXC_WEIGHT[0]

# MAX_EXC_WEIGHT = np.arange(20+de, 30, de)

GNA_mean = []
GNA_var = []
GNA_fano = []

parameters2D = {
    "L": 3.0, # Length of the square domain (mm) or diameter of the circular domain (mm)
    "rho": 71, # Neuron density (neurons/mm^2)
    "soma_diameter": 0.015, # Diameter of the soma (mm)
    "d_mean": 0.15, # Dendritic tree diameter (mm) Gaussian distribution
    "d_sigma": 0.02, # Standard deviation of the dendritic tree diameter (mm)
    "l_mean": 1.0, # Mean axon length (mm) Rayleigh distribution
    "segment_length": 0.01, # Segment length for axon growth (mm)
    "sigma_axon_angle": 0.1, # Standard deviation of the axon angle (radians)
    "agg": 1, # Control parameter for the aggregation of neurons
    "n_centers": 10, # Number of centers for the aggregation of neurons
    "base_sigma": 0.8, # Base standard deviation for aggregation centers (mm)
    "BC_type": "SW",  # Boundary conditions type (PBC: Periodic Boundary Conditions, SW: Sticky Walls)
    "geometry": "square", # Geometry of the domain (square or circle)
    "seed": 0 # Random seed for reproducibility
}

alpha = 0.2

if parameters2D["geometry"] == "circle" and parameters2D["BC_type"] == "PBC":
    raise ValueError("Periodic Boundary Conditions (PBC) are not compatible with a circular geometry. Please choose 'SW' for BC_type or change the geometry to 'square'.")
if parameters2D["geometry"] == "square":
    N = int(parameters2D["rho"] * parameters2D["L"]**2) # Total number of neurons
elif parameters2D["geometry"] == "circle":
    R = parameters2D["L"] / np.sqrt(np.pi)
    print(f"Effective radius of the circular domain: {R:.2f} mm")
    N = int(parameters2D["rho"] * np.pi * R**2) # Total number of neurons, same as in the square case for the same L, but we keep it here for clarity

# print(f"Number of neurons: {N}")
parametersDynamics = {
    "max_exc_weight": 18,
    "max_inh_weight": 0.5,
    "noise_max": 3,
    "sim_time": 10000,
    "seed": 42
}

new_sim = True
sims = 10

cell_size = 0.2
n_grid = int(parameters2D["L"] / cell_size)

for sim in range(sims):

    output_file = f"GNA_map_vsGini/GNAmap_{parameters2D['geometry']}_{parameters2D['BC_type']}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_sim{sim}_ch{chunk_id}.txt"
    open(output_file, "w").close()  # Limpiar el archivo antes de escribir
    with open(output_file, "w") as f:
        f.write("base_sigma\tMAX_EXC_WEIGHT\tGini\tmean_peaks_height\tstd_peaks_height\tfano_factor\n")

    for s in base_sigmas:

        parameters2D["base_sigma"] = s
        filename_neuron_params = fcn.run_2D_tune_positions(parameters2D, new_sim=new_sim)
        filename_network_tries = fcn.create_network2D(parameters2D, new_sim=new_sim)
        Gini = fcn.measure_gini2D(filename_neuron_params, parameters2D, n_grid)

        A, filename_adj = fcn.modify_A_alpha2D(
            fcn.load_A(filename_network_tries),
            alpha,
            parameters2D,
            new_sim=new_sim,
            sim=sim
        )

        N = A.shape[0]

        for e in MAX_EXC_WEIGHT:

            parametersDynamics["max_exc_weight"] = e
            parametersDynamics["out_filename"] = filename_adj

            filename_dynamics = fcn.run_dynamics_2D(
                parametersDynamics,
                parameters2D,
                alpha,
                sim_number=sim,
                new_sim=new_sim, 
            )

            firings_t, firings_i = fcn.load_firings(filename_dynamics)

            # print("N usado en Python:", N)
            # print("Máxima neurona en firings:", firings_i.max())

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

            # inter_func_corr = fcn.compute_interfunctional_community_GNA_correlation(
            #     filename_dynamics=filename_dynamics,
            #     SIM_TIME=parametersDynamics["sim_time"],
            #     N=N,
            #     W=25,
            #     min_community_size=5
            # )

            GNA_mean.append(GNA.mean())
            GNA_var.append(GNA.var())
            GNA_fano.append(fano_factor)

            with open(output_file, "a") as f:
                f.write(f"{s:.4f}\t{e:.2f}\t{Gini:.6f}\t{mean_peaks_height:.6f}\t{std_peaks_height:.6f}\t{fano_factor:.6f}\n")
            
            # print(f"base_sigma: {s:.4f}, MAX_EXC_WEIGHT: {e:.2f}, Gini: {Gini:.6f}, mean_peaks_height: {mean_peaks_height:.6f}, std_peaks_height: {std_peaks_height:.6f}, fano_factor: {fano_factor:.6f}")
        print(f"Completada base_sigma: {s:.4f} para sim {sim}")

# Le cuesta 10h
# Mirar distribuciones de kIn con las dos condiciones periodicas de contorno
# Plotear las neuronas y cambiar e tamaño del nodo en funcion del grado para cada condicones, hacer smoothing de los colores
# Poner en el mapa ruster en extreos e intemedios
# Grupo de neuronas inactivas
# Hacer el cultivo cilindrico con SW en las paredes y en las bases que se corte el crecimiento. 
# Mirar el small_world, richclub  