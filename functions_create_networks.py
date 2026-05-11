from fileinput import filename
import subprocess
import numpy as np
from pathlib import Path
from pypalettes import load_cmap
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
import networkx as nx
from numba import njit

sunset2 = load_cmap('Sunset2', cmap_type='continuous')

# @njit
def load_neuron_params(filename, nrows_skiped=15):
        data = np.loadtxt(filename, skiprows=nrows_skiped)
        return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
# # @njit
def run_2D_tune_positions(parameters2D, show_prints_c = False, show_prints_control = True, new_sim = False):

    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['geometry']}_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['geometry']}_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"

    filepath = Path(filename)

    # =========================
    # COMPROBAR SI EXISTE
    # =========================

    if filepath.exists() and (new_sim == False):
        if show_prints_control:
            print("Configuration already exists:", filename)

    else:
        if show_prints_control:
            print("Creating configuration:", filename)

        result = subprocess.run([
            "./2D_tune_positions.exe",
            str(parameters2D["L"]),
            str(parameters2D["rho"]),
            str(parameters2D["soma_diameter"]),
            str(parameters2D["d_mean"]),
            str(parameters2D["d_sigma"]),
            str(parameters2D["l_mean"]),
            str(parameters2D["segment_length"]),
            str(parameters2D["sigma_axon_angle"]),
            str(parameters2D["agg"]),
            str(parameters2D["n_centers"]),
            str(parameters2D["base_sigma"]),
            str(parameters2D["BC_type"]),
            str(parameters2D["geometry"])
        ], capture_output=True, text=True)

        if show_prints_c:
            print(result.stdout)
            print(result.stderr)

    return filename
# @njit
def create_network2D(parameters2D, show_prints_c = False, show_prints_control = True, new_sim = False):
    
    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['geometry']}_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['geometry']}_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"

    # --------------------------------

    if parameters2D["agg"] == 0:
        output_filename = f"2D_created_networks/2D_adjacency_matrix_{parameters2D['geometry']}_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        output_filename = f"2D_created_networks/2D_adjacency_matrix_{parameters2D['geometry']}_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"

    filepath = Path(output_filename)
    if filepath.exists() and (new_sim == False):
        if show_prints_control:
            print("Configuration already exists:", output_filename)
    else:
        if show_prints_control:
            print("Creating configuration:", output_filename)

        result = subprocess.run(["./2D_create_network.exe", filename],
        capture_output=True,
        text=True)
        
        if show_prints_c:
            print(result.stdout)
            print(result.stderr)

    return output_filename
# @njit
def run_dynamics_2D(parametersDynamics, parameters2D, alpha, show_prints_c=False, show_prints_control=True, sim_number=0, new_sim=False):

    if parameters2D["agg"] == 0:
        filename = f"2D_dynamics/2D_dynamics_{parameters2D['geometry']}_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}_sim{sim_number}.txt"
    else:
        filename = f"2D_dynamics/2D_dynamics_{parameters2D['geometry']}_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}_sim{sim_number}.txt"

    path = Path(filename)
    if path.exists() and (new_sim == False):
        if show_prints_control:
            print("Dynamics already exists:", filename)
    else:
        if show_prints_control:
            print("Simulating dynamics:", filename)

        result = subprocess.run([
            "./2D_dynamics.exe",
            parametersDynamics["out_filename"],
            str(parametersDynamics["max_exc_weight"]),
            str(parametersDynamics["max_inh_weight"]),
            str(parametersDynamics["noise_max"]),
            str(parametersDynamics["sim_time"]),
            str(parametersDynamics["seed"]),
            str(sim_number),
            parameters2D["geometry"]
        ], capture_output=False, text=False)

        if show_prints_c:
            print(result.stdout)
            print(result.stderr)
    
    return filename

# @njit
def run_dynamics_3D(parametersDynamics, parameters3D, alpha, show_prints_c=False, show_prints_control=True, new_sim=False, sim_number=0):

    if parameters3D["agg"] == 0:
        filename = f"3D_dynamics/3D_dynamics_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}_sim{sim_number}.txt"
    else:
        filename = f"3D_dynamics/3D_dynamics_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}_sim{sim_number}.txt"

    path = Path(filename)
    if path.exists() and (new_sim == False):
        if show_prints_control:
            print("Dynamics already exists:", filename)
    else:
        if show_prints_control:
            print("Simulating dynamics:", filename)

        result = subprocess.run([
            "./3D_dynamics.exe",
            parametersDynamics["out_filename"],
            str(parametersDynamics["max_exc_weight"]),
            str(parametersDynamics["max_inh_weight"]),
            str(parametersDynamics["noise_max"]),
            str(parametersDynamics["sim_time"]),
            str(parametersDynamics["seed"]),
            str(sim_number)
        ], capture_output=True, text=True)

        if show_prints_c:
            print(result.stdout)
            print(result.stderr)
    
    return filename
    
# @njit
def load_firings(path):
    # Cargar ignorando la primera línea (# time neuron)
    data = np.loadtxt(path, comments="#")

    # Separar columnas
    firings_t = data[:, 0].astype(int)
    firings_i = data[:, 1].astype(int)

    return firings_t, firings_i
# @njit
def neuron_map_gini2D(filename, parameters2D):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')
    viridis = load_cmap('Viridis', cmap_type='continuous')

    X, Y, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params(filename)

    L = parameters2D["L"]
    xmin = ymin = -L/2
    xmax = ymax = L/2
    n_grid = 20

    # =========================
    # HISTOGRAMA 2D
    # =========================

    H, _, _ = np.histogram2d(X, Y, bins=n_grid, range=[[xmin,xmax],[ymin,ymax]])

    # =========================
    # GINI
    # =========================

    bin_sum = H.sum()

    sum_gini = 0.0
    H_flat = H.flatten()
    H_flat.sort()
    H_flat = H_flat[::-1]

    G = np.zeros(len(H_flat))

    for i, h in enumerate(H_flat):
        sum_gini += h
        G[i] = sum_gini/bin_sum

    bins = np.arange(1, n_grid*n_grid+1)
    bins = bins/(n_grid*n_grid)

    area_between = np.trapezoid(G-bins, bins)
    Gini = 2*area_between

    print("Gini coefficient:", Gini)

    # =========================
    # FIGURA CON DOS SUBPLOTS
    # =========================

    fig, axs = plt.subplots(1, 2, figsize=(12,5))

    # -------------------------
    # RED DE NEURONAS
    # -------------------------

    ax = axs[0]


    ax.set_xlabel('X (mm)', fontsize=16, labelpad=10)
    ax.set_ylabel('Y (mm)', fontsize=16, labelpad=10)

    ax.set_xlim(-L/2, L/2)
    ax.set_ylim(-L/2, L/2)

    # ax.set_xticks([])
    # ax.set_yticks([])

    x_grid = np.linspace(xmin, xmax, n_grid+1)
    y_grid = np.linspace(ymin, ymax, n_grid+1)

    for x in x_grid:
        ax.axvline(x, color='black', lw=0.3)

    for y in y_grid:
        ax.axhline(y, color='black', lw=0.3)

    ax.imshow(
        H.T,
        extent=[xmin,xmax,ymin,ymax],
        origin='lower',
        aspect='equal',
        alpha=0.5,
        cmap='cividis',
    )

    ax.scatter(X, Y, s=7, color=viridis(0.05), zorder=3)

    ax.tick_params(axis='both', which='major', labelsize = 14)
    ax.tick_params(axis='both', which='minor', labelsize = 0)

    ax.set_title("Neuron positions")

    # -------------------------
    # CURVA DE GINI
    # -------------------------

    ax= axs[1]

    ax.plot(bins, G, color=sunset2(0.1), marker='o', markersize=3, lw=1)
    ax.plot(bins, bins, color='black', lw=2)

    ax.legend(fontsize=15, loc='lower right', frameon=False)
    # ax.text(0.75, 0.15, r"$\Lambda = $" + f"{Gini:.3f}", fontsize=30, ha='center', va='center', transform=ax.transAxes)

    ax.set_ylabel('Occupation', fontsize=18, labelpad=10)

    ax.tick_params(axis='both', which='major', labelsize = 14)
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.text(0.58, 0.5, 'Perfect equality', fontsize = 18, color='black', ha='center', va='center', rotation=45, transform=ax.transAxes)


    ax.set_xlabel('Cumulative grid cells', fontsize=18, labelpad=10)

    # ax.set_title(f"Gini coefficient = {Gini:.3f}", fontsize=16)

    plt.tight_layout()
    plt.savefig("plots/2D_neuron_map_gini.svg", dpi=300, bbox_inches='tight')
    plt.show()
# @njit
def run_3D_tune_positions(parameters3D, show_prints_c = False, show_prints_control = True, new_sim = False):
    if parameters3D["agg"] == 0:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    else:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    filepath = Path(filename)

    # =========================
    # COMPROBAR SI EXISTE
    # =========================

    if filepath.exists() and not new_sim:
        if show_prints_control:
            print("Configuration already exists:", filename)

    else:
        if show_prints_control:
            print("Creating configuration:", filename)

        result = subprocess.run([
            "./3D_tune_positions.exe",
            str(parameters3D["L"]),
            str(parameters3D["rho"]),
            str(parameters3D["soma_diameter"]),
            str(parameters3D["d_mean"]),
            str(parameters3D["d_sigma"]),
            str(parameters3D["l_mean"]),
            str(parameters3D["segment_length"]),
            str(parameters3D["sigma_pol"]),
            str(parameters3D["sigma_azi"]),
            str(parameters3D["agg"]),
            str(parameters3D["n_centers"]),
            str(parameters3D["base_sigma"]),
            str(parameters3D["BC_type"])
        ], capture_output=True, text=True)

        if show_prints_c:
            print(result.stdout)
            print(result.stderr)
        
    return filename
# @njit
def neuron_map_adj2D(filename_in, A, parameters2D, n_grid=20):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

    # =========================
    # CARGA DE DATOS
    # =========================
    X, Y, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params(filename_in)

    # =========================
    # DOMINIO ESPACIAL
    # =========================
    L = parameters2D["L"]
    xmin = ymin = -L / 2
    xmax = ymax = L / 2

    H, _, _ = np.histogram2d(
        X, Y,
        bins=n_grid,
        range=[[xmin, xmax], [ymin, ymax]]
    )

    # =========================
    # FIGURA
    # =========================
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # -------------------------
    # IZQUIERDA: MAPA ESPACIAL
    # -------------------------
    ax = axs[0]

    x_grid = np.linspace(xmin, xmax, n_grid + 1)
    y_grid = np.linspace(ymin, ymax, n_grid + 1)

    for x in x_grid:
        ax.axvline(x, color='black', lw=0.3, zorder=1)
    for y in y_grid:
        ax.axhline(y, color='black', lw=0.3, zorder=1)

    ax.imshow(
        H.T,
        extent=[xmin, xmax, ymin, ymax],
        origin='lower',
        aspect='equal',
        alpha=0.5,
        cmap=sunset2,
        zorder=0
    )

    ax.scatter(X, Y, s=7, color='firebrick', zorder=3)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('X (mm)', fontsize=16, labelpad=10)
    ax.set_ylabel('Y (mm)', fontsize=16, labelpad=10)
    ax.tick_params(axis='both', which='major', labelsize=14)

    # -------------------------
    # DERECHA: MATRIZ DE ADYACENCIA
    # -------------------------
    ax = axs[1]

    ax.imshow(A, aspect='equal', interpolation='nearest', cmap='viridis')
    ax.set_xlim(0, A.shape[0])
    ax.set_ylim(A.shape[0], 0)
    ax.set_xlabel('neuron', fontsize=16, labelpad=10)
    ax.set_ylabel('neuron', fontsize=16, labelpad=10)
    ax.tick_params(axis='both', which='major', labelsize=14)

    plt.tight_layout()
    plt.show()
# @njit
def create_network3D(parameters3D, show_prints_c = False, show_prints_control = True, new_sim = False):

    if parameters3D["agg"] == 0:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    else:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"

    # --------------------------------

    if parameters3D["agg"] == 0:
        output_filename = f"3D_created_networks/3D_adjacency_matrix_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    else:
        output_filename = f"3D_created_networks/3D_adjacency_matrix_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    
    filepath = Path(output_filename)
    if filepath.exists() and not new_sim:
        if show_prints_control:
            print("Configuration already exists:", output_filename)
    else:
        if show_prints_control:
            print("Creating configuration:", output_filename)

        result = subprocess.run(["./3D_create_network.exe", filename],
        capture_output=True,
        text=True)
        
        if show_prints_c:
            print(result.stdout)
            print(result.stderr)

    return output_filename
# @njit
def load_neuron_params_3D(filename, nrows_skiped=15):
    data = np.loadtxt(filename, skiprows=nrows_skiped)
    return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
# @njit
def neuron_map_gini3D(filename, parameters3D):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

    # filename = '3D_initial_configurations/3D_neurons_params_agg_nc50_s0.10_L2.0_rho125_l1.00.txt'

    # ====== CARGA ======
    X, Y, Z, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params_3D(filename)

    # ====== DOMINIO ======
    L = parameters3D["L"]
    xmin = ymin = zmin = -L/2
    xmax = ymax = zmax =  L/2

    # ====== HISTOGRAMA 3D ======
    n_grid = 20

    H, edges = np.histogramdd(
        sample=np.column_stack([X, Y, Z]),
        bins=(n_grid, n_grid, n_grid),
        range=((xmin, xmax), (ymin, ymax), (zmin, zmax))
    )

    bin_sum = H.sum()

    H_flat = H.flatten().astype(np.float64)
    H_flat.sort()
    H_flat = H_flat[::-1]

    G = np.zeros(H_flat.size)

    sum_gini = 0.0
    for i, h in enumerate(H_flat):
        sum_gini += h
        G[i] = sum_gini/bin_sum

    bins = np.arange(1, H_flat.size + 1)
    bins = bins/H_flat.size

    area_between = np.trapezoid(G-bins, bins)
    Gini = 2*area_between

    # print("Gini coefficient (3D):", Gini)

    # =========================
    # FIGURA CON DOS SUBPLOTS
    # =========================

    fig = plt.figure(figsize=(12,6))
    # gs = fig.add_gridspec(1, 2, width_ratios=[2, 4])

    # --------- RED 3D ---------

    ax1 = fig.add_subplot(121, projection='3d')

    ax1.scatter(
        X, Y, Z,
        s=8,
        color='firebrick',
        depthshade=False
    )

    ax1.set_xlim(-L/2, L/2)
    ax1.set_ylim(-L/2, L/2)
    ax1.set_zlim(-L/2, L/2)

    ax1.set_xlabel('X (mm)', fontsize=16, labelpad=10)
    ax1.set_ylabel('Y (mm)', fontsize=16, labelpad=10)
    ax1.set_zlabel('Z (mm)', fontsize=16, labelpad=10)

    ax1.tick_params(labelsize=14)

    ax1.set_box_aspect([1,1,1])

    # ax1.set_title("Neuron positions (3D)")

    # --------- GINI ---------

    ax2 = fig.add_subplot(122)

    ax2.plot(
        bins,
        G,
        color=sunset2(0.1),
        marker='o',
        markersize=2.5,
        linestyle='-',
        lw=1
    )

    ax2.plot(
        bins,
        bins,
        color='black',
        lw=2,
        label="Perfect equality"
    )

    ax2.legend(fontsize=14, loc='lower right', frameon=False)
    ax2.text(0.75, 0.15, r"$\Delta = $" + f"{Gini:.3f}", fontsize=20, ha='center', va='center')

    ax2.set_ylabel('Occupation', fontsize=16, labelpad = 10)

    ax2.tick_params(labelsize=14)

    # ax2.set_title(f"Gini coefficient = {Gini:.3f}")

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    plt.show()       


#### DYAMICS NETWORK ANALYSIS ####

# @njit
def load_A(path):
    A = np.loadtxt(path, dtype=np.uint8)
    np.fill_diagonal(A, 0)
    return A
# @njit
def connectivity_from_A(A: np.ndarray):
    A = np.asarray(A)
    assert A.ndim == 2 and A.shape[0] == A.shape[1], "A debe ser cuadrada"
    N = A.shape[0]

    # por si acaso
    A = A.copy()
    np.fill_diagonal(A, 0)

    E = A.sum()                 # número de enlaces (dirigidos)
    C = E / (N * (N - 1))       # densidad
    k_mean = E / N              # grado medio (out/in)

    return C, k_mean, E

import numpy as np
import matplotlib.pyplot as plt
from pypalettes import load_cmap
# @njit
def load_neuron_params_3D(filename, nrows_skiped=15):
    data = np.loadtxt(filename, skiprows=nrows_skiped)
    return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
# @njit
def load_A(path):
    A = np.loadtxt(path, dtype=np.uint8)
    np.fill_diagonal(A, 0)
    return A
# @njit
def neuron_map_adj3D(filename_in, A, parameters3D):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

    # =========================
    # CARGA
    # =========================
    X, Y, Z, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params_3D(filename_in)
    # A = load_A(filename_out)

    # =========================
    # DOMINIO
    # =========================
    L = parameters3D["L"]

    # =========================
    # FIGURA
    # =========================
    fig = plt.figure(figsize=(12, 6))

    # --------- IZQUIERDA: RED 3D ---------
    ax1 = fig.add_subplot(121, projection='3d')

    ax1.scatter(
        X, Y, Z,
        s=8,
        color='firebrick',
        depthshade=False
    )

    ax1.set_xlim(-L/2, L/2)
    ax1.set_ylim(-L/2, L/2)
    ax1.set_zlim(-L/2, L/2)

    ax1.set_xlabel('X (mm)', fontsize=16, labelpad=10)
    ax1.set_ylabel('Y (mm)', fontsize=16, labelpad=10)
    ax1.set_zlabel('Z (mm)', fontsize=16, labelpad=10)

    ax1.tick_params(labelsize=14)
    ax1.set_box_aspect([1, 1, 1])

    # fija la vista para que no quede raro
    ax1.view_init(elev=20, azim=-60)

    # --------- DERECHA: MATRIZ DE ADYACENCIA ---------
    ax2 = fig.add_subplot(122)

    ax2.imshow(
        A,
        aspect='equal',
        interpolation='nearest',
    )

    ax2.set_xlim(0, A.shape[0])
    ax2.set_ylim(A.shape[0], 0)

    ax2.set_xlabel('neuron', fontsize=16, labelpad=10)
    ax2.set_ylabel('neuron', fontsize=16, labelpad=10)

    ax2.tick_params(labelsize=14)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    plt.show()
from scipy.signal import find_peaks

# =========================================================
# FUNCIÓN PARA CALCULAR GNA A PARTIR DEL RASTER
# =========================================================
def compute_GNA_from_raster(firings_t, firings_i, N, SIM_TIME, W=25):
    
    firings_t = np.asarray(firings_t, dtype=int)
    firings_i = np.asarray(firings_i, dtype=int) - 1  # 0-based
    
    # Para cada neurona, guardamos el último tiempo en que disparó
    last_spike = -np.inf * np.ones(N)
    
    GNA = np.zeros(SIM_TIME)
    active_count = 0
    
    # Ordenar por tiempo (importante)
    order = np.argsort(firings_t)
    firings_t = firings_t[order]
    firings_i = firings_i[order]
    
    idx = 0
    
    for t in range(1, SIM_TIME + 1):
        
        # Añadir spikes en este tiempo
        while idx < len(firings_t) and firings_t[idx] == t:
            neuron = firings_i[idx]
            
            # Si estaba inactiva, ahora cuenta
            if last_spike[neuron] < t - W:
                active_count += 1
            
            last_spike[neuron] = t
            idx += 1
        
        # Quitar neuronas que han dejado de estar activas
        active_count = np.sum(last_spike >= t - W)
        
        GNA[t-1] = active_count / N
    
    return GNA
# @njit
def find_GNA_peaks(GNA):
    prom = 0.2 * np.std(GNA)

    peaks, properties = find_peaks(
        GNA,
        prominence=prom,  # prominencia mínima del pico
        height = 0.05,
        distance=20
    )

    peak_heights = properties.get("peak_heights", GNA[peaks] if len(peaks) > 0 else np.array([]))

    if len(peak_heights) == 0:
        mean_peak_height = 0.0
        std_peak_height = 0.0
    else:
        mean_peak_height = peak_heights.mean()
        std_peak_height = peak_heights.std()
    
    return peaks, peak_heights, mean_peak_height, std_peak_height
# def compute_GNA_from_raster(firings_t, firings_i, N, SIM_TIME, W=25, indices_are_1_based=True):
#     """
#     Calcula la actividad global de red (GNA) como:
#         GNA(t) = nº de neuronas activas en la ventana / N

#     donde una neurona cuenta solo UNA vez por ventana,
#     aunque dispare varios spikes dentro de ella.

#     Parámetros
#     ----------
#     firings_t : array
#         Tiempos de disparo (en ms, enteros)
#     firings_i : array
#         Índices de neurona asociados a cada spike
#     N : int
#         Número total de neuronas
#     SIM_TIME : int
#         Tiempo total de simulación (ms)
#     W : int
#         Tamaño de la ventana deslizante (ms)
#     indices_are_1_based : bool
#         True si firings_i va de 1..N
#         False si firings_i va de 0..N-1

#     Devuelve
#     --------
#     GNA : array de longitud SIM_TIME
#         Fracción de neuronas activas en cada instante
#     GNA_percent : array de longitud SIM_TIME
#         Igual que GNA pero en %
#     active_neurons_count : array de longitud SIM_TIME
#         Número de neuronas activas en cada instante
#     """

#     firings_t = np.asarray(firings_t, dtype=int)
#     firings_i = np.asarray(firings_i, dtype=int)

#     # Pasar índices a 0-based si están en 1..N
#     if indices_are_1_based:
#         neuron_ids = firings_i - 1
#     else:
#         neuron_ids = firings_i.copy()

#     # Filtrar por seguridad spikes fuera de rango
#     valid = (
#         (firings_t >= 1) & (firings_t <= SIM_TIME) &
#         (neuron_ids >= 0) & (neuron_ids < N)
#     )
#     firings_t = firings_t[valid]
#     neuron_ids = neuron_ids[valid]

#     # -----------------------------------------------------
#     # 1) Matriz binaria de actividad: activity[n, t] = 1 si
#     #    la neurona n dispara en el instante t
#     # -----------------------------------------------------
#     activity = np.zeros((N, SIM_TIME), dtype=np.uint8)

#     # firings_t va de 1..SIM_TIME, así que restamos 1 para indexar 0..SIM_TIME-1
#     activity[neuron_ids, firings_t - 1] = 1

#     # -----------------------------------------------------
#     # 2) Para cada neurona, calculamos si ha estado activa
#     #    al menos una vez dentro de la ventana de tamaño W
#     #
#     #    Usamos convolución temporal por neurona:
#     #    si la suma en la ventana > 0, esa neurona está activa
#     # -----------------------------------------------------
#     kernel = np.ones(W, dtype=int)

#     active_in_window = np.zeros_like(activity, dtype=np.uint8)

#     for n in range(N):
#         counts_in_window = np.convolve(activity[n], kernel, mode='same')
#         active_in_window[n] = (counts_in_window > 0).astype(np.uint8)

#     # -----------------------------------------------------
#     # 3) Número de neuronas activas por instante
#     # -----------------------------------------------------
#     active_neurons_count = np.sum(active_in_window, axis=0)

#     # Fracción y porcentaje
#     GNA = active_neurons_count / N
#     # GNA_percent = 100 * GNA

#     return GNA, active_neurons_count

# @njit
def Gini_coefficient(X, Y, n_grid, xmin, xmax, ymin, ymax):
    H, _, _ = np.histogram2d(X, Y, bins=n_grid, range=[[xmin,xmax],[ymin,ymax]])
    bin_sum = H.sum()

    sum_gini = 0.0
    H_flat = H.flatten()
    H_flat.sort()
    H_flat = H_flat[::-1]

    G = np.zeros(len(H_flat))

    for i, h in enumerate(H_flat):
        sum_gini += h
        G[i] = sum_gini/bin_sum

    bins = np.arange(1, len(H_flat)+1)
    bins = bins/len(H_flat)

    area_between = np.trapezoid(G-bins, bins)
    Gini = 2*area_between
    return Gini

# @njit
def modify_A_alpha2D(A, alpha, parameters2D, seed=None, new_sim = False):
    """
    A: matriz de conteos (número de intersecciones m_ij)
    alpha: probabilidad de éxito por intento
    parameters2D: diccionario con los parámetros de la red 2D (necesarios para calcular la probabilidad efectiva)
    seed: semilla para reproducibilidad (opcional)
    new_sim: flag para indicar si es una nueva simulación
    """
    if parameters2D["agg"] == 0:
        output_filename = f"2D_dynamics/networks_for_dynamics/2D_adjacency_matrix_{parameters2D['geometry']}_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}.txt"
    else:
        output_filename = f"2D_dynamics/networks_for_dynamics/2D_adjacency_matrix_{parameters2D['geometry']}_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}.txt"

    path = Path(output_filename)
    if path.exists() and not new_sim:
        # print(f"Modified adjacency matrix already exists with alpha = {alpha:.2f}:", output_filename)
        A_new = load_A(output_filename)
        A_new = get_GC(A_new)
    else:
        if not (0 <= alpha <= 1):
            raise ValueError("alpha debe estar entre 0 y 1")

        rng = np.random.default_rng(seed)

        # Probabilidad efectiva por par (i,j)
        P = 1.0 - (1.0 - alpha) ** A

        # Sorteo Bernoulli con esa probabilidad
        A_new = (rng.random(A.shape) < P).astype(np.uint8)

        # Quitar autoconexiones
        np.fill_diagonal(A_new, 0)
        A_new = get_GC(A_new)

        pd.DataFrame(A_new).to_csv(output_filename, sep=' ', index=False, header=False)
    

    return A_new, output_filename


# @njit
def modify_A_alpha3D(A, alpha, parameters3D, seed=None, new_sim = False):
    """
    A: matriz de conteos (número de intersecciones m_ij)
    alpha: probabilidad de éxito por intento
    parameters3D: diccionario con los parámetros de la red 3D (necesarios para calcular la probabilidad efectiva)
    seed: semilla para reproducibilidad (opcional)
    new_sim: flag para indicar si es una nueva simulación
    """
    if parameters3D["agg"] == 0:
        output_filename = f"3D_dynamics/networks_for_dynamics/3D_adjacency_matrix_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}_a{alpha:.2f}.txt"
    else:
        output_filename = f"3D_dynamics/networks_for_dynamics/3D_adjacency_matrix_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}_a{alpha:.2f}.txt"

    path = Path(output_filename)
    if path.exists() and not new_sim:
        print(f"Modified adjacency matrix already exists with alpha = {alpha:.2f}:", output_filename)
        A_new = load_A(output_filename)
    else:
        if not (0 <= alpha <= 1):
            raise ValueError("alpha debe estar entre 0 y 1")

        rng = np.random.default_rng(seed)

        # Probabilidad efectiva por par (i,j)
        P = 1.0 - (1.0 - alpha) ** A

        # Sorteo Bernoulli con esa probabilidad
        A_new = (rng.random(A.shape) < P).astype(np.uint8)

        # Quitar autoconexiones
        np.fill_diagonal(A_new, 0)
        A_new = get_GC(A_new)

        pd.DataFrame(A_new).to_csv(output_filename, sep=' ', index=False, header=False)
    
    return A_new, output_filename
# @njit
def detect_peaks_adaptive(signal, time=None, smooth=True,
                          window_length=11, polyorder=3,
                          prominence_factor=0.25,
                          min_distance=None,
                          min_width=None):
    """
    Detecta picos de forma adaptativa usando prominencia relativa
    a la escala de la propia señal.

    Parámetros
    ----------
    signal : array-like
        Señal de entrada.
    time : array-like or None
        Vector de tiempos. Si se proporciona, min_distance y min_width
        se interpretan en las mismas unidades que time.
    smooth : bool
        Si True, aplica suavizado Savitzky-Golay antes de detectar.
    window_length : int
        Tamaño de ventana del suavizado (debe acabar siendo impar).
    polyorder : int
        Orden del polinomio para Savitzky-Golay.
    prominence_factor : float
        Factor multiplicativo sobre una escala robusta de la señal
        (p90 - p10) para fijar la prominencia mínima.
    min_distance : float or None
        Distancia mínima entre picos, en unidades de time si time se da.
    min_width : float or None
        Anchura mínima de pico, en unidades de time si time se da.

    Devuelve
    --------
    peaks : np.ndarray
        Índices de los picos detectados.
    properties : dict
        Propiedades devueltas por scipy.signal.find_peaks.
    x_used : np.ndarray
        Señal realmente usada para detectar (suavizada o no).
    """

    x = np.asarray(signal, dtype=float).copy()

    if x.ndim != 1:
        raise ValueError("signal debe ser un array 1D")

    if len(x) < 3:
        return np.array([], dtype=int), {}, x

    # 1. Suavizado opcional
    if smooth:
        wl = min(window_length, len(x))
        if wl % 2 == 0:
            wl -= 1
        if wl < 3:
            wl = 3

        po = min(polyorder, wl - 1)

        if wl >= 3 and po >= 1:
            x = savgol_filter(x, window_length=wl, polyorder=po)

    median = np.median(x)
    mad = np.median(np.abs(x - median))
    scale = mad

    # Si la señal es casi constante
    if mad <= 0:
        return np.array([], dtype=int), {}, x

    # 3. Conversión de min_distance y min_width a muestras
    distance_samples = 1
    width_samples = None

    if time is not None:
        time = np.asarray(time, dtype=float)
        if len(time) != len(x):
            raise ValueError("time y signal deben tener la misma longitud")

        dt = time[1] - time[0]
        if dt <= 0:
            raise ValueError("time debe ser estrictamente creciente")

        if min_distance is not None:
            distance_samples = max(1, int(min_distance / dt))

        if min_width is not None:
            width_samples = max(1, int(min_width / dt))
    else:
        if min_distance is not None:
            distance_samples = max(1, int(min_distance))
        if min_width is not None:
            width_samples = max(1, int(min_width))

    baseline = np.percentile(x, 30)

    peaks, properties = find_peaks(
        x,
        prominence=5 * mad,
        height=baseline + 4 * mad,
        distance=distance_samples,
        width=width_samples
    )

    return peaks, properties, x

# @njit
def build_temporal_neuron_series(filename, time_window=25, SIM_TIME = 5000, N=1000):

    data = np.loadtxt(filename, skiprows=1)
    firings_t = data[:, 0].astype(int)
    firings_i = data[:, 1].astype(int)


    neurons_ids = np.arange(1, N+1)  # IDs de neuronas de 1 a N
    times = np.arange(0, SIM_TIME, time_window)

    # Quiero pasar de (t,i) a una matriz M de tamaño (len(times), N) donde M[j,n] = 1 si la neurona n ha disparado al menos una vez en el intervalo [times[j], times[j]+time_window)    
    M = np.zeros((len(times), N), dtype=np.uint8)

    for t, neuron_id in zip(firings_t, firings_i):

        row = t // time_window
        col = neuron_id - 1  # neuronas de 1..N a índices de 0..N-1

        if neuron_id < 1 or neuron_id > N:
            continue
        if row >= len(times):
            continue

        M[row, col] = 1

    return M, times
# @njit
def pearson_correlation(x, y):

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("x e y deben tener la misma forma")

    if x.ndim != 1:
        raise ValueError("x e y deben ser vectores 1D")

    n = len(x)
    if n < 2:
        return np.nan

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    numerator = np.sum((x - mean_x) * (y - mean_y))
    denominator = np.sqrt(np.sum((x - mean_x) ** 2) * np.sum((y - mean_y) ** 2))

    if denominator == 0:
        return np.nan

    r = numerator / denominator
    return r

# @njit
def pearson_corr_matrix(M):

    N = M.shape[1]
    corr_matrix = np.zeros((N, N), dtype=float)

    for i in range(N):
        for j in range(i, N):
            r = pearson_correlation(M[:, i], M[:, j])
            corr_matrix[i, j] = r
            corr_matrix[j, i] = r

    return corr_matrix
# @njit
def plot_raster_interleaved_inhibitory(
    firings_t,
    firings_i,
    ax=None,
    inhibitory_fraction=0.20,
    color=None,
    s=7,
    xlabel='Time (ms)',
    ylabel='Neuron',
    fontsize=16,
    ticksize=14
):
    """
    Plotea un raster reorganizando las neuronas inhibitorias para que no aparezcan
    todas agrupadas arriba, sino intercaladas entre las excitadoras.

    Supone que:
    - los índices de neurona empiezan en 1
    - las últimas inhibitory_fraction*N neuronas son inhibitorias

    Parámetros
    ----------
    firings_t : array-like
        Tiempos de disparo.
    firings_i : array-like
        Índices de neurona (1-based).
    ax : matplotlib.axes.Axes, opcional
        Eje sobre el que dibujar. Si es None, crea figura nueva.
    inhibitory_fraction : float
        Fracción de neuronas inhibitorias. Por defecto 0.20.
    color : color matplotlib, opcional
        Color de los puntos.
    s : float
        Tamaño de los puntos del scatter.
    """

    firings_t = np.asarray(firings_t)
    firings_i = np.asarray(firings_i, dtype=int)

    if firings_i.ndim != 1 or firings_t.ndim != 1:
        raise ValueError("firings_t y firings_i deben ser arrays 1D.")
    if len(firings_t) != len(firings_i):
        raise ValueError("firings_t y firings_i deben tener la misma longitud.")
    if np.any(firings_i < 1):
        raise ValueError("Se esperan índices de neurona empezando en 1.")

    # Número total de neuronas inferido del raster
    N = int(np.max(firings_i))

    # Separación excitadoras / inhibitorias
    N_exc = int(round((1.0 - inhibitory_fraction) * N))
    N_inh = N - N_exc

    # Índices originales en base 0
    exc = np.arange(0, N_exc)          # excitadoras
    inh = np.arange(N_exc, N)          # inhibitorias

    # Construimos un nuevo orden intercalando inhibitorias uniformemente
    new_order = []
    ratio = N_exc / N_inh if N_inh > 0 else np.inf

    e_idx = 0
    i_idx = 0

    while e_idx < N_exc or i_idx < N_inh:
        target_exc = int(round((i_idx + 1) * ratio)) if N_inh > 0 else N_exc

        while e_idx < min(target_exc, N_exc):
            new_order.append(exc[e_idx])
            e_idx += 1

        if i_idx < N_inh:
            new_order.append(inh[i_idx])
            i_idx += 1

    while e_idx < N_exc:
        new_order.append(exc[e_idx])
        e_idx += 1

    new_order = np.array(new_order, dtype=int)

    # Mapa: índice viejo -> índice nuevo
    old_to_new = np.empty(N, dtype=int)
    old_to_new[new_order] = np.arange(N)

    # Reindexar firings_i (1-based -> 0-based -> reordenado -> 1-based)
    firings_i_reordered = old_to_new[firings_i - 1] + 1

    return firings_i_reordered, new_order

from numba import njit
# @njit
def degree_distribution(A_binary, normalize=False):
    degrees = np.sum(A_binary, axis=1).astype(int)
    k_vals = np.arange(degrees.max() + 1)

    pk = np.array([
        np.sum(degrees == k) for k in k_vals
    ], dtype=float)

    if normalize and pk.sum() > 0:
        pk /= pk.sum()

    return k_vals, pk

# @njit
def binarize_functional_matrix(F, percentile=90):
    """
    Binariza una matriz funcional F utilizando un umbral basado en percentiles.

    Procedimiento:
    1. Se ignora la diagonal (autoconexiones) asignando NaN.
    2. Se consideran solo los valores de la parte triangular superior
       (para evitar duplicados en matrices simétricas).
    3. Se calcula un umbral como el percentil indicado (por defecto, 90).
       → El percentil p es el valor por debajo del cual se encuentra el p% de los datos.
    4. Se construye una matriz binaria donde:
       - 1 si F_ij >= threshold (conexiones más fuertes)
       - 0 en caso contrario
    5. Se simetriza la matriz para obtener una red no dirigida.
    6. Se eliminan autoconexiones (diagonal = 0).

    Interpretación:
    - Con percentile=90, se conservan aproximadamente el 10% de las conexiones más fuertes.
    - El umbral es relativo a la distribución de F, por lo que controla la densidad
      de la red de forma adaptativa.

    Parámetros
    ----------
    F : ndarray
        Matriz funcional (por ejemplo, correlaciones).
    percentile : float
        Percentil usado para definir el umbral (0–100).

    Returns
    -------
    F_binary : ndarray
        Matriz binaria (0/1) de conectividad.
    threshold : float
        Valor del umbral aplicado.
    """
    F_tmp = F.copy()
    np.fill_diagonal(F_tmp, np.nan)

    iu = np.triu_indices_from(F_tmp, k=1)
    vals = F_tmp[iu]
    valid_vals = vals[~np.isnan(vals)]
    F_binary = np.zeros_like(F, dtype=np.uint8)


    # threshold = np.percentile(valid_vals, percentile)
    # F_binary[iu] = (F[iu] >= threshold).astype(np.uint8)

    threshold = np.percentile(np.abs(valid_vals), percentile)
    F_binary[iu] = (np.abs(F[iu]) >= threshold).astype(np.uint8)

    F_binary = F_binary + F_binary.T
    np.fill_diagonal(F_binary, 0)

    return F_binary, threshold

# @njit
def order_F_louvain_from_raster(
    F,
    F_binary,
    firings_i,
    seed=0
):
    rng = np.random.default_rng(seed)
    N = F.shape[0]

    active = np.unique(firings_i.astype(int) - 1)
    all_neurons = np.arange(N)
    inactive = np.setdiff1d(all_neurons, active)

    print(f"Neuronas activas: {len(active)}")
    print(f"Neuronas inactivas: {len(inactive)}")

    W_active = F_binary[np.ix_(active, active)].copy()
    np.fill_diagonal(W_active, 0)

    G = nx.from_numpy_array(W_active)

    communities = nx.community.louvain_communities(
        G,
        weight="weight",
        resolution=1.0,
        seed=seed
    )

    communities = sorted(communities, key=len, reverse=True)

    degree = np.sum(W_active, axis=1)

    order_active_local = []

    for comm in communities:
        comm = list(comm)

        comm_sorted = sorted(
            comm,
            key=lambda i: degree[i],
            reverse=True
        )

        order_active_local.extend(comm_sorted)

    order_active_local = np.array(order_active_local)
    order_active = active[order_active_local]

    final_order = np.empty(N, dtype=int)

    inactive_positions = rng.choice(
        N,
        size=len(inactive),
        replace=False
    )

    mask_inactive = np.zeros(N, dtype=bool)
    mask_inactive[inactive_positions] = True

    final_order[mask_inactive] = rng.permutation(inactive)
    final_order[~mask_inactive] = order_active

    F_ordered = F[np.ix_(final_order, final_order)]
    np.fill_diagonal(F_ordered, 0)

    return F_ordered, final_order, communities


# @njit
def rewire_directed_fast_inplace(A, n_rewires):
    N = A.shape[0]

    for _ in range(n_rewires):

        while True:
            i = np.random.randint(0, N)
            j = np.random.randint(0, N)

            if i != j and A[i, j] == 1:
                break

        while True:
            u = np.random.randint(0, N)
            v = np.random.randint(0, N)

            if u != v and A[u, v] == 0:
                break

        A[i, j] = 0
        A[u, v] = 1

# @njit
def plot_full_panel(
    F_ordered,
    A_ordered,
    firings_t,
    firings_i_louvain,
    GNA,
    k_vals_A,
    pk_A,
    k_vals,
    pk,
    k_vals_rw,
    pk_rw,
    SIM_TIME,
    base_sigma,
    MAX_EXC_WEIGHT,
    savepath=None,
    N = 2000
):
    fig = plt.figure(figsize=(14, 6))

    gs = fig.add_gridspec(
        2, 3,
        width_ratios=[1, 1, 2.8],
        height_ratios=[1, 1],
        wspace=0.28,
        hspace=0.28
    )

    # ========================================================
    # A) Functional matrix - arriba izquierda
    # ========================================================

    axA = fig.add_subplot(gs[0, 0])

    axA.imshow(
        F_ordered,
        aspect="equal",
        interpolation="nearest",
        cmap="binary",
        vmin=0,
        vmax=1
    )

    axA.set_ylim(F_ordered.shape[0], 0)
    axA.set_title("A) Functional", fontsize=14)
    axA.set_xticks([0, 1000, 2000])
    axA.set_xticklabels(["0", "1k", "2k"])
    axA.set_yticks([0, 1000, 2000])
    axA.set_yticklabels(["0", "1k", "2k"])
    axA.tick_params(labelsize=10)

    # ========================================================
    # B) Structural matrix - abajo izquierda
    # ========================================================

    axB = fig.add_subplot(gs[1, 0])

    axB.imshow(
        A_ordered,
        aspect="equal",
        interpolation="nearest",
        cmap="binary",
        vmin=0,
        vmax=1
    )

    axB.set_ylim(A_ordered.shape[0], 0)
    axB.set_title("B) Structural", fontsize=14)
    axB.set_xticks([0, 1000, 2000])
    axB.set_xticklabels(["0", "1k", "2k"])
    axB.set_yticks([0, 1000, 2000])
    axB.set_yticklabels(["0", "1k", "2k"])
    axB.tick_params(labelsize=10)

    # ========================================================
    # C) Raster - arriba centro
    # ========================================================

    axC = fig.add_subplot(gs[0, 1])

    axC.scatter(
        firings_t / 1000,
        firings_i_louvain,
        s=0.2,
        color=sunset2(0.1)
    )

    axC.set_ylim(N, 0)

    axC.set_title("C) Raster", fontsize=14)
    axC.set_ylabel("Neuron", fontsize=12)
    axC.set_xlim(0, 2)
    axC.set_yticks([0, 1000, 2000])
    axC.set_yticklabels(["0", "1k", "2k"])
    axC.tick_params(labelsize=10)

    # ========================================================
    # D) GNA - abajo centro
    # ========================================================

    axD = fig.add_subplot(gs[1, 1])

    axD.plot(
        np.arange(1, SIM_TIME + 1) / 1000,
        GNA,
        lw=1.5,
        color=sunset2(0.1)
    )

    axD.set_title("D) GNA", fontsize=14)
    axD.set_xlabel("Time (s)", fontsize=12)
    axD.set_ylabel("GNA", fontsize=12)
    axD.set_xlim(0, 2)
    # axD.set_ylim(0, 0.78)
    axD.tick_params(labelsize=10)

    # ========================================================
    # E) Degree distribution - derecha grande
    # ========================================================

    axE = fig.add_subplot(gs[:, 2])

    axE.plot(
        k_vals,
        pk,
        lw=2,
        color=sunset2(0.1),
        label="Functional"
    )

    axE.plot(
        k_vals_rw,
        pk_rw,
        lw=2,
        color=sunset2(0.5),
        label="Rewired"
    )

    axE.plot(
        k_vals_A,
        pk_A,
        lw=2,
        color=sunset2(0.9),
        label="Structural"
    )

    axE.set_title("E) Degree distribution", fontsize=14)
    axE.set_xlabel("Degree (k)", fontsize=12)
    axE.set_ylabel("P(k)", fontsize=12)
    axE.legend(frameon=False, fontsize=10)
    axE.tick_params(labelsize=10)
    axE.set_xlim(0, 400)
    fig.suptitle(
        rf"$\sigma={base_sigma:.2f}$, $E_{{exc}}={MAX_EXC_WEIGHT:.2f}$",
        fontsize=16,
        y=0.98
    )

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")

    plt.show()


def compile_c_codes():
    gcc = "C:/msys64/ucrt64/bin/gcc.exe"

    subprocess.run([
        gcc, 
        "-O3",
        "2D_create_network.c", 
        "2D_functions.c",
        "-o", 
        "2D_create_network.exe",
        "-lm"
    ], check=True)

    subprocess.run([
        gcc, 
        "-Wall", 
        "-O3",
        "2D_tune_positions.c", 
        "2D_functions.c",
        "-o", 
        "2D_tune_positions.exe",
        "-lm"
    ], check=True)

    subprocess.run([
        gcc,
        "-O3",
        "2D_dynamics.c", 
        "2D_functions.c",
        "-o", 
        "2D_dynamics.exe",
        "-lm"
    ], check=True)
# @njit
def get_GC(A):
    G = nx.from_numpy_array(A, create_using=nx.DiGraph)

    largest_cc = max(nx.strongly_connected_components(G), key=len)

    GC = G.subgraph(largest_cc).copy()

    A_GC = nx.to_numpy_array(GC, dtype=np.uint8)

    return A_GC

def create_network(parametersStructure, alpha = 0.2, seed=None, new_sim=False, build_network = True, show_prints_c = False, show_prints_control = False):
    
    filename_neuron_params = run_2D_tune_positions(parametersStructure, new_sim=new_sim, show_prints_c=show_prints_c, show_prints_control=show_prints_control)
    
    filename_adj = None
    
    if build_network:
        ADJ_matrix, filename_adj = modify_A_alpha2D(
            load_A(create_network2D(parametersStructure, new_sim=new_sim, show_prints_c=show_prints_c, show_prints_control=show_prints_control)),
            alpha,
            parametersStructure,
            seed=seed,
            new_sim=new_sim
        )

    return filename_neuron_params, filename_adj

def create_network_3D(parametersStructure, alpha = 0.2, seed=None, new_sim=False, build_network = True, show_prints_c = False, show_prints_control = False):
    
    filename_neuron_params = run_3D_tune_positions(parametersStructure, new_sim=new_sim, show_prints_c=show_prints_c, show_prints_control=show_prints_control)
    
    filename_adj = None
    
    if build_network:
        ADJ_matrix, filename_adj = modify_A_alpha3D(
            load_A(create_network3D(parametersStructure, new_sim=new_sim, show_prints_c=show_prints_c, show_prints_control=show_prints_control)),
            alpha,
            parametersStructure,
            seed=seed,
            new_sim=new_sim
        )

    return filename_neuron_params, filename_adj

# @njit
def plot_degree_distribution(k_vals_F, pk_F, k_vals_rw, pk_rw, k_vals_A, pk_A, savepath=None):
    plt.figure(figsize=(6, 5))

    mean_degree_F = np.sum(k_vals_F * pk_F)
    plt.plot(
        k_vals_F,
        pk_F,
        lw=2,
        color=viridis(0.1),
        label="Functional"+r"$\langle k \rangle = {:.0f}$".format(mean_degree_F)
    )

    plt.plot(
        k_vals_rw,
        pk_rw,
        lw=2,
        color=viridis(0.5),
        label="Rewired Functional"
    )

    mean_degree_A = np.sum(k_vals_A * pk_A)
    plt.plot(
        k_vals_A,
        pk_A,
        lw=2,
        color=viridis(0.9),
        label="Structural"+r"$\langle k \rangle = {:.0f}$".format(mean_degree_A)
    )

    # plt.title("Degree distribution", fontsize=14)
    plt.xlabel("Degree (k)", fontsize=18)
    plt.ylabel("P(k)", fontsize=18)
    plt.legend(frameon=False, fontsize=14)
    plt.tick_params(labelsize=16)
    # plt.xlim(0, 400)

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")

    plt.show()
# @njit
def build_F_matrix(filename_dynamics, SIM_TIME=5000, N=1000):
    M, times = build_temporal_neuron_series(
        filename_dynamics,
        time_window=25,
        SIM_TIME=SIM_TIME,
        N=N
    )

    firings_t, firings_i = np.loadtxt(filename_dynamics, skiprows=1, unpack=True)

    F = np.corrcoef(M, rowvar=False)

    F_binary, threshold = binarize_functional_matrix(F, percentile=90)

    F_ordered, order, communities = order_F_louvain_from_raster(
        F,
        F_binary,
        firings_i,
        seed=0
    )

    return F, F_binary, F_ordered, order, communities

# @njit
def rewire_and_compute_degree_distribution(A_binary, n_rewires, seed=None):
    # Number of rewires is equal to the number of edges in the original direted network
    n_rewires = n_rewires if n_rewires is not None else np.sum(A_binary)
    A_rewired = A_binary.copy()
    rewire_directed_fast_inplace(A_rewired, n_rewires)
    k_vals_rw, pk_rw = degree_distribution(A_rewired, normalize=True)
    return A_rewired, k_vals_rw, pk_rw

# @njit
def properties_of_network(A, directed=True):

    """
    Calculates properties of the network such as:
    -  Degree distribution
    -   Clustering coefficient
    -   Average path length
    -   Centrality measures
    -   Modularity
    -   Mean degree
    -   Degree assortativity
    -   Etc.

    Parameters
    ----------
    A : ndarray
        Adjacency matrix of the network.
    Returns
    -------
    properties : dict
        Dictionary containing the calculated properties.
  
    """

    G = nx.from_numpy_array(A, create_using=nx.DiGraph if directed else nx.Graph)

    properties = {}

    properties['mean_degree'] = np.mean(np.sum(A, axis=1))
    properties['degree_distribution'] = degree_distribution(A, normalize=True)
    properties['clustering_coefficient'] = nx.average_clustering(G.to_undirected())
    properties['average_path_length'] = nx.average_shortest_path_length(G) if nx.is_connected(G.to_undirected()) else np.inf
    properties['degree_assortativity'] = nx.degree_assortativity_coefficient(G)
    properties['modularity'] = None  # Requires community detection

    return properties

viridis = load_cmap('Viridis', cmap_type='continuous')
# @njit
def plot_spatial_distr_AND_adjacency(filename_neuron_params, filename_adj, parametersStructure, order=None):
    
    X, Y, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params(filename_neuron_params)

    L = parametersStructure["L"]
    geometry = parametersStructure["geometry"]

    xmin = ymin = -L/2
    xmax = ymax = L/2
    n_grid = 30

    R = L / 2

    H, x_edges, y_edges = np.histogram2d(
        X, Y,
        bins=n_grid,
        range=[[xmin, xmax], [ymin, ymax]]
    )

    # =====================================================
    # Máscara geométrica para el cálculo del Gini
    # =====================================================

    if geometry == "square":
        H_valid = H.flatten()

    elif geometry == "circle":
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

        Xc, Yc = np.meshgrid(x_centers, y_centers, indexing="ij")

        mask_circle = Xc**2 + Yc**2 <= R**2

        H_valid = H[mask_circle]

    else:
        raise ValueError(f"Unknown geometry: {geometry}")

    # =====================================================
    # Gini
    # =====================================================

    bin_sum = H_valid.sum()

    H_flat = np.sort(H_valid.flatten())[::-1]

    G = np.zeros(len(H_flat))
    sum_gini = 0.0

    for i, h in enumerate(H_flat):
        sum_gini += h
        G[i] = sum_gini / bin_sum

    bins = np.arange(1, len(H_flat) + 1)
    bins = bins / len(H_flat)

    area_between = np.trapezoid(G - bins, bins)
    Gini = 2 * area_between

    # =====================================================
    # PLOTS
    # =====================================================

    fig, axs = plt.subplots(
        1, 2,
        figsize=(12, 5),
        gridspec_kw={"width_ratios": [1, 1], "wspace": 0.3}
    )

    # =====================================================
    # Neurons in space
    # =====================================================

    ax1 = axs[0]

    ax1.set_xlabel('X (mm)', fontsize=18, labelpad=10)
    ax1.set_ylabel('Y (mm)', fontsize=18, labelpad=10)

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.set_aspect("equal")

    x_grid = np.linspace(xmin, xmax, n_grid + 1)
    y_grid = np.linspace(ymin, ymax, n_grid + 1)

    for x in x_grid:
        ax1.axvline(x, color='black', lw=0.3, alpha=0.4)

    for y in y_grid:
        ax1.axhline(y, color='black', lw=0.3, alpha=0.4)

    H_plot = H.T.copy()

    if geometry == "circle":
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        Xc, Yc = np.meshgrid(x_centers, y_centers, indexing="xy")

        mask_circle = Xc**2 + Yc**2 <= R**2

        H_plot = np.ma.masked_where(~mask_circle, H_plot)

    ax1.imshow(
        H_plot,
        extent=[xmin, xmax, ymin, ymax],
        origin='lower',
        aspect='equal',
        alpha=0.5,
        cmap='cividis',
    )

    ax1.scatter(X, Y, s=7, color=viridis(0.05), zorder=3)

    if geometry == "circle":
        circle = plt.Circle(
            (0, 0),
            R,
            fill=False,
            color="black",
            lw=1.5,
            zorder=4
        )
        ax1.add_patch(circle)
    
    ax1.text(
        0.05, 0.95,
        r"$\Lambda = $" + f"{Gini:.2f}",
        transform=ax1.transAxes,
        fontsize=14,
        verticalalignment='top',
        bbox=dict(
            boxstyle="square,pad=0.55",
            facecolor="white",
            alpha=0.9,
            edgecolor="none"
        ), zorder = 5
        )
        
    ax1.tick_params(axis='both', which='major', labelsize=16)

    # =====================================================
    # Adjacency matrix
    # =====================================================

    A = load_A(filename_adj)

    if order is not None:
        A = A[np.ix_(order, order)]

    ax2 = axs[1]

    ax2.imshow(
        A,
        aspect="equal",
        interpolation="nearest",
        cmap="binary",
        vmin=0,
        vmax=1
    )

    ax2.set_xlabel('Neuron', fontsize=18)
    ax2.set_ylabel('Neuron', fontsize=18)
    ax2.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()
# @njit
def compute_spatial_gini(X, Y, L, geometry="square", n_grid=30):
    """
    Calcula un índice tipo Gini de heterogeneidad espacial usando una cuadrícula.

    Para geometry='square':
        usa todos los bins del dominio cuadrado [-L/2, L/2] x [-L/2, L/2].

    Para geometry='circle':
        usa solo los bins cuyo centro cae dentro del círculo de radio R = L/2.

    Returns
    -------
    Gini : float
        Índice de heterogeneidad espacial.
    H : np.ndarray
        Histograma 2D completo.
    x_edges, y_edges : np.ndarray
        Bordes de los bins.
    mask_valid : np.ndarray
        Máscara de bins válidos para el cálculo del Gini.
    """

    xmin = ymin = -L / 2
    xmax = ymax =  L / 2

    H, x_edges, y_edges = np.histogram2d(
        X, Y,
        bins=n_grid,
        range=[[xmin, xmax], [ymin, ymax]]
    )

    if geometry == "square":

        mask_valid = np.ones_like(H, dtype=bool)

    elif geometry == "circle":

        R = L / 2

        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

        Xc, Yc = np.meshgrid(x_centers, y_centers, indexing="ij")

        mask_valid = Xc**2 + Yc**2 <= R**2

    else:
        raise ValueError(f"Unknown geometry: {geometry}")

    H_valid = H[mask_valid]

    if H_valid.sum() == 0:
        return np.nan, H, x_edges, y_edges, mask_valid

    H_flat = np.sort(H_valid.flatten())[::-1]

    cumulative = np.cumsum(H_flat) / H_flat.sum()

    bins = np.arange(1, len(H_flat) + 1)
    bins = bins / len(H_flat)

    area_between = np.trapezoid(cumulative - bins, bins)

    Gini = 2 * area_between

    return Gini, H, x_edges, y_edges, mask_valid
# @njit
def plot_F_rusterGNA(filename_dynamics, F, order, SIM_TIME = 10000, N=1000, xlim = None):

    firings_t, firings_i = load_firings(filename_dynamics)

    GNA = compute_GNA_from_raster(firings_t, firings_i, N = N, SIM_TIME=SIM_TIME, W = 25)

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

    F_ordered = F[np.ix_(order, order)]

    fig = plt.figure(figsize=(8, 5))

    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1., 1],
        height_ratios=[3,1],
        wspace=0.40,
        hspace=0.02
    )

    # =====================================================
    # Izquierda: matriz funcional
    # =====================================================
    ax_matrix = fig.add_subplot(gs[:, 0])

    ax_matrix.imshow(
        F_ordered,
        aspect="equal",
        interpolation="nearest",
        cmap="binary",
        vmin=0,
        vmax=1
    )

    ax_matrix.set_ylim(F_ordered.shape[0], 0)
    ax_matrix.set_xlabel("Neurons", fontsize=18)
    ax_matrix.set_ylabel("Neurons", fontsize=18)
    ax_matrix.tick_params(axis='both', which='major',
                          labelsize=14, length=6, width=1.2)

    # =====================================================
    # Derecha arriba: raster plot
    # =====================================================

    #Ordenar el raster con el mismo orden que la matriz funcional
    old_to_new = np.empty(N, dtype=int)
    old_to_new[order] = np.arange(N)

    firings_i = old_to_new[firings_i - 1] + 1

    ax_raster = fig.add_subplot(gs[0, 1])

    ax_raster.scatter(
        firings_t,
        firings_i,
        s=0.2,
        color=sunset2(0.1)
    )

    ax_raster.set_ylabel("Neurons", fontsize=18)
    ax_raster.tick_params(axis='both', which='major',
                          labelsize=14, length=6, width=1.2)
    ax_raster.tick_params(axis='x', labelbottom=False)

    ax_raster.set_ylim(N, 0)
    ax_raster.set_xlim(0, xlim if xlim is not None else SIM_TIME)

    # =====================================================
    # Derecha abajo: GNA
    # =====================================================
    ax_gna = fig.add_subplot(gs[1, 1], sharex=ax_raster)

    ax_gna.plot(
        np.arange(1, SIM_TIME + 1),
        GNA,
        lw=1.5,
        color=sunset2(0.1)
    )

    ax_gna.set_xlabel("Time (s)", fontsize=18)
    ax_gna.set_ylabel("GNA", fontsize=18)
    ax_gna.tick_params(axis='both', which='major',
                       labelsize=14, length=6, width=1.2)
    
    ax_gna.set_xlim(0, xlim if xlim is not None else SIM_TIME)

    # Opcional: limitar ventana temporal
    # ax_gna.set_xlim(0, 2)

    plt.tight_layout()
    plt.show()

# Función para detectar nodos centrales usando centralidad de grado
# @njit
def detect_central_nodes(A, top_k=10):
    G = nx.from_numpy_array(A, create_using=nx.DiGraph)
    degree_centrality = nx.degree_centrality(G)
    sorted_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
    central_nodes = sorted_nodes[:top_k]
    return central_nodes


def compute_network_metrics_from_GC(A_GC, N_total, seed=0):
    """
    Calcula métricas sobre la componente gigante.

    A_GC : np.ndarray
        Matriz de adyacencia de la componente gigante.
    N_total : int
        Número total de neuronas antes de extraer la GC.
    """

    n_gc = A_GC.shape[0]

    if n_gc <= 1:
        return {
            "n_gc": n_gc,
            "gc_fraction": n_gc / N_total,
            "mean_degree": np.nan,
            "density": np.nan,
            "modularity": np.nan,
            "clustering": np.nan,
            "avg_shortest_path": np.nan,
            "global_efficiency": np.nan,
            "diameter": np.nan,
        }

    G = nx.from_numpy_array(A_GC, create_using=nx.DiGraph)
    G_und = G.to_undirected()

    mean_degree_in = np.mean([d for _, d in G.in_degree()])
    mean_degree_out = np.mean([d for _, d in G.out_degree()])
    density = nx.density(G)

    gc_fraction = n_gc / N_total

    clustering = nx.average_clustering(G_und)

    communities = nx.community.louvain_communities(
        G_und,
        seed=seed
    )

    modularity = nx.community.modularity(
        G_und,
        communities
    )

    avg_shortest_path = nx.average_shortest_path_length(G_und)

    global_efficiency = nx.global_efficiency(G_und)

    diameter = nx.diameter(G_und)

    assortativity = nx.degree_assortativity_coefficient(G_und)

    return {
        "n_gc": n_gc,
        "gc_fraction": gc_fraction,
        "mean_degree_in": mean_degree_in,
        "mean_degree_out": mean_degree_out,
        "density": density,
        "modularity": modularity,
        "clustering": clustering,
        "avg_shortest_path": avg_shortest_path,
        "global_efficiency": global_efficiency,
        "diameter": diameter,
        "assortativity": assortativity
    }

