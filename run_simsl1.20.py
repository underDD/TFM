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

n_chunks = 5
chunk_id = int(sys.argv[1])

# fcn.compile_c_codes3D()

base_sigmas = [0.04]
axonal_lengths = np.linspace(0.1, 1.5, 50)
axonal_lengths = np.array_split(axonal_lengths, n_chunks)[chunk_id]

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
    "l_mean": 0.3, # Mean axon length (mm) Rayleigh distribution
    "segment_length": 0.01, # Segment length for axon growth (mm)
    "sigma_pol": 0.1, # Standard deviation of the axon angle (radians)
    "sigma_azi": 0.1, # Standard deviation of the azimuthal
    "agg": 1, # Control parameter for the aggregation of neurons
    "n_centers": 400, # Number of centers for the aggregation of neurons
    "base_sigma": 0.8, # Base standard deviation for aggregation centers (mm)
    "BC_type": "SW",  # Boundary conditions type (PBC: Periodic Boundary Conditions, SW: Sticky Walls)
    "geometry": "cube", # Geometry of the domain (cube or cylinder)
    "seed": 12345
}

alpha = 0.2

if parameters3D["geometry"] == "cylinder" and parameters3D["BC_type"] == "PBC":
    raise ValueError("Periodic Boundary Conditions (PBC) are not compatible with a cylindrical geometry. Please choose 'SW' for BC_type or change the geometry to 'cube'.")
if parameters3D["geometry"] == "cube":
    N = int(parameters3D["rho"] * parameters3D["L"]**3) # Total number of neurons
elif parameters3D["geometry"] == "cylinder":
    R = parameters3D["L"] / 2
    print(f"Effective radius of the cylindrical domain: {R:.2f} mm")
    N = int(parameters3D["rho"] * np.pi * R**2) # Total number of neurons, same as in the square case for the same L, but we keep it here for clarity

# print(f"Number of neurons: {N}")
parametersDynamics = {
    "max_exc_weight": 12,
    "max_inh_weight": 0.5,
    "noise_max": 3,
    "sim_time": 10000,
    "seed": 42
}

new_sim = True
sims = 1
# I want to tho 6 simulations more becaouse i have don e already 4 starting from 0
for sim in range(sims):

    for s in base_sigmas:
        parameters3D["base_sigma"] = s
        parameters3D["seed"] = 12345 + sim*1000  # Update seed for each simulation
        filename_neuron_params = fcn.run_3D_tune_positions(parameters3D, new_sim=False)
        
        output_file = f"Gini_things/3D_GNAmap_{parameters3D['geometry']}_{parameters3D['BC_type']}_rho{parameters3D['rho']}_s{parameters3D['base_sigma']:.2f}_sim{sim}_chunk{chunk_id}.txt"
        open(output_file, "w").close()  # Limpiar el archivo antes de escribir
        with open(output_file, "w") as f:
            f.write("axonal_length\tMAX_EXC_WEIGHT\tmean_peaks_height\tstd_peaks_height\tfano_factor\tGC_fraction\tmean_degree\n")
        
        for l in axonal_lengths:

            parameters3D["l_mean"] = l
            filename_network_tries = fcn.create_network3D(parameters3D, new_sim=new_sim)

            A, filename_adj = fcn.modify_A_alpha3D(
                fcn.load_A(filename_network_tries),
                alpha,
                parameters3D,
                seed = 12345 + sim*1000,
                new_sim=new_sim,
                sim=sim
            )

            N_GC = A.shape[0]
            GC_fraction = N_GC / N
            mean_degree_in = A.sum(axis=0).mean()

            parametersDynamics["out_filename"] = filename_adj

            for e in MAX_EXC_WEIGHT:

                parametersDynamics["max_exc_weight"] = e

                filename_dynamics = fcn.run_dynamics_3D(
                    parametersDynamics,
                    parameters3D,
                    alpha,
                    sim_number=sim,
                    new_sim=new_sim, 
                    show_prints_control=False
                )

                firings_t, firings_i = fcn.load_firings(filename_dynamics)

                GNA = fcn.compute_GNA_from_raster(
                        firings_t=firings_t,
                        firings_i=firings_i,
                        N=N_GC,
                        SIM_TIME=parametersDynamics["sim_time"],
                        W=25
                    )
                peaks, peak_heights, mean_peaks_height, std_peaks_height = fcn.find_GNA_peaks(GNA)
                fano_factor = std_peaks_height**2 / mean_peaks_height if mean_peaks_height > 0 else 0.0

                GNA_mean.append(GNA.mean())
                GNA_var.append(GNA.var())
                GNA_fano.append(fano_factor)

                with open(output_file, "a") as f:
                    f.write(f"{l:.2f}\t{e:.2f}\t{mean_peaks_height:.6f}\t{std_peaks_height:.6f}\t{fano_factor:.6f}\t{GC_fraction:.6f}\t{mean_degree_in:.2f}\n")
                
                print(f"base_sigma: {s:.2f}, axonal_length: {l:.2f}, MAX_EXC_WEIGHT: {e:.2f}, mean_peaks_height: {mean_peaks_height:.6f}, std_peaks_height: {std_peaks_height:.6f}, fano_factor: {fano_factor:.6f}, GC_fraction: {GC_fraction:.6f}, mean_degree_in: {mean_degree_in:.2f}")

    print(f"Finished simulation {sim+1}/{sims}")