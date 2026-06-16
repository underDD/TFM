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

import inspect

print(fcn.__file__)
print(inspect.getsource(fcn.load_neuron_params_3D))

# fcn.compile_c_codes3D()

calibration_filename = "figures_TFM/data_for_figures/gini3D_n_grid13_calibration.txt"
data = np.loadtxt(calibration_filename, skiprows=1)
gini_calibration = data[:, 0]
base_sigmas_calibration = data[:, 1]

n_chunks = 4
chunk_id = int(sys.argv[1])

base_sigmas = np.array_split(base_sigmas_calibration, n_chunks)[chunk_id]

MAX_EXC_WEIGHT = np.linspace(2, 20.0, 50)
# de = MAX_EXC_WEIGHT[1] - MAX_EXC_WEIGHT[0]

# MAX_EXC_WEIGHT = np.arange(20+de, 30, de)

GNA_mean = []
GNA_var = []
GNA_fano = []

parameters3D = {
    "L": 2.0, # Length of the square domain (mm) or diameter of the circular domain (mm)
    "height": 1.0, # Height of the cylinder domain (mm), only used if geometry is "cylinder"
    "rho": 1400, # Neuron density (neurons/mm^2)
    "soma_diameter": 0.015, # Diameter of the soma (mm)
    "d_mean": 0.15, # Dendritic tree diameter (mm) Gaussian distribution
    "d_sigma": 0.02, # Standard deviation of the dendritic tree diameter (mm)
    "l_mean": 0.7, # Mean axon length (mm) Rayleigh distribution
    "segment_length": 0.01, # Segment length for axon growth (mm)
    "sigma_pol": 0.1, # Standard deviation of the axon angle (radians)
    "sigma_azi": 0.1, # Standard deviation of the azimuthal angle (radians)
    "agg": 1, # Control parameter for the aggregation of neurons
    "n_centers": 112, # Number of centers for the aggregation of neurons
    "base_sigma": 0.04, # Base standard deviation for aggregation centers (mm)
    "BC_type": "SW",  # Boundary conditions type (PBC: Periodic Boundary Conditions, SW: Sticky Walls)
    "geometry": "cubeLongRange", # Geometry of the domain (cube or sphere)
    "seed": 0,
}

if parameters3D["geometry"] == "cylinder" and parameters3D["BC_type"] == "PBC":
    raise ValueError("Periodic Boundary Conditions (PBC) are not compatible with a cylindrical geometry. Please choose 'SW' for BC_type or change the geometry to 'cube'.")
if parameters3D["geometry"] == "cube" or parameters3D["geometry"] == "cubeLongRange" or parameters3D["geometry"] == "cubeBiased":
    N = int(parameters3D["rho"] * parameters3D["L"]**2 * parameters3D["height"]) # Total number of neurons
elif parameters3D["geometry"] == "cylinder":
    R = parameters3D["L"] / np.sqrt(np.pi)
    print(f"Effective radius of the cylinder domain: {R:.2f} mm")
    N = int(parameters3D["rho"] * np.pi * R**2 * parameters3D["height"]) # Total number of neurons, same as in the cube case for the same L, but we keep it here for clarity

alpha = 0.2
# print(f"Number of neurons: {N}")
parametersDynamics = {
    "max_exc_weight": 8,
    "max_inh_weight": 0.5,
    "noise_max": 3,
    "sim_time": 10000,
    "seed": 42
}

new_sim = True
sims = 10

cell_size = 0.15
n_grid = int(parameters3D["L"] / cell_size)

for sim in range(sims):

    output_file = f"GNA_map_vsGini/3D_GNAmap_{parameters3D['geometry']}_{parameters3D['BC_type']}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_sim{sim}_ch{chunk_id}.txt"
    open(output_file, "w").close()  # Limpiar el archivo antes de escribir
    with open(output_file, "w") as f:
        f.write("base_sigma\tMAX_EXC_WEIGHT\tGini\tmean_peaks_height\tstd_peaks_height\tfano_factor\n")

    for s in base_sigmas:

        parameters3D["base_sigma"] = s
        filename_neuron_params = fcn.run_3D_tune_positions(parameters3D, new_sim=new_sim)
        filename_network_tries = fcn.create_network3D(parameters3D, new_sim=new_sim)
        X,Y, Z, _,_ = fcn.load_neuron_params_3D(filename_neuron_params)
        Gini = fcn.compute_gini_3D(X, Y, Z, parameters3D["L"], parameters3D["height"], geometry=parameters3D["geometry"], cell_size = cell_size)
        M = fcn.load_A3D(filename_network_tries, N)
        A, filename_adj, nodes_GSCC = fcn.modify_A_alpha3D_fast(M, alpha, parameters3D=parameters3D, new_sim=True)


        N_A = A.shape[0]

        for e in MAX_EXC_WEIGHT:

            parametersDynamics["max_exc_weight"] = e
            parametersDynamics["out_filename"] = filename_adj

            filename_dynamics = fcn.run_dynamics_3D(
                parametersDynamics,
                parameters3D,
                alpha,
                new_sim=new_sim, 
                show_prints_c = False, 
                show_prints_control = False
            )

            firings_t, firings_i = fcn.load_firings(filename_dynamics)

            # print("N usado en Python:", N)
            # print("Máxima neurona en firings:", firings_i.max())

            assert firings_i.max() <= N_A

            GNA = fcn.compute_GNA_from_raster(
                    firings_t=firings_t,
                    firings_i=firings_i,
                    N=N_A,
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
        print(f"Completada base_sigma: {s:.4f} para sim {sim}")

# Le cuesta 10h
# Mirar distribuciones de kIn con las dos condiciones periodicas de contorno
# Plotear las neuronas y cambiar e tamaño del nodo en funcion del grado para cada condicones, hacer smoothing de los colores
# Poner en el mapa ruster en extreos e intemedios
# Grupo de neuronas inactivas
# Hacer el cultivo cilindrico con SW en las paredes y en las bases que se corte el crecimiento. 
# Mirar el small_world, richclub  