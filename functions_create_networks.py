import subprocess
import numpy as np
from pathlib import Path
from pypalettes import load_cmap
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks, savgol_filter


def load_neuron_params(filename, nrows_skiped=14):
        data = np.loadtxt(filename, skiprows=nrows_skiped)
        return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

def run_2D_tune_positions(parameters2D, show_prints_c = False, show_prints_control = True):

    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"

    filepath = Path(filename)

    # =========================
    # COMPROBAR SI EXISTE
    # =========================

    if filepath.exists():
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
            str(parameters2D["BC_type"])
        ], capture_output=True, text=True)

        if show_prints_c:
            print(result.stdout)
            print(result.stderr)

    return filename

def create_network2D(parameters2D, show_prints_c = True, show_prints_control = True):
    
    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"

    # --------------------------------

    if parameters2D["agg"] == 0:
        output_filename = f"2D_created_networks/2D_adjacency_matrix_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    else:
        output_filename = f"2D_created_networks/2D_adjacency_matrix_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}.txt"
    
    filepath = Path(output_filename)
    if filepath.exists():
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

def run_dynamics_2D(parametersDynamics, parameters2D, alpha, show_prints=True):

    if parameters2D["agg"] == 0:
        filename = f"2D_dynamics/2D_dynamics_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}.txt"
    else:
        filename = f"2D_dynamics/2D_dynamics_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}_exc{parametersDynamics['max_exc_weight']:.3f}.txt"

    path = Path(filename)
    if path.exists():
        if show_prints:
            print("Dynamics already exists:", filename)
    else:
        result = subprocess.run([
            "./2D_dynamics.exe",
            parametersDynamics["out_filename"],
            str(parametersDynamics["max_exc_weight"]),
            str(parametersDynamics["max_inh_weight"]),
            str(parametersDynamics["noise_max"]),
            str(parametersDynamics["sim_time"]),
            str(parametersDynamics["seed"])
        ], capture_output=True, text=True)

        if show_prints:
            print(result.stdout)
            print(result.stderr)
    
    return filename
    
def load_firings(path):
    # Cargar ignorando la primera línea (# time neuron)
    data = np.loadtxt(path, comments="#")

    # Separar columnas
    firings_t = data[:, 0].astype(int)
    firings_i = data[:, 1].astype(int)

    return firings_t, firings_i

def neuron_map_gini2D(filename, parameters2D):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

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

    fig, axs = plt.subplots(1, 2, figsize=(12,6))

    # -------------------------
    # RED DE NEURONAS
    # -------------------------

    ax = axs[0]


    ax.set_xlabel('X (mm)', fontsize=16, labelpad=10)
    ax.set_ylabel('Y (mm)', fontsize=16, labelpad=10)

    ax.set_xlim(-L/2, L/2)
    ax.set_ylim(-L/2, L/2)

    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())

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
        cmap=sunset2
    )

    ax.scatter(X, Y, s=7, color='firebrick', zorder=3)
    ax.tick_params(axis='both', which='major', labelsize = 14)

    # ax.set_title("Neuron positions")

    # -------------------------
    # CURVA DE GINI
    # -------------------------

    ax = axs[1]

    ax.plot(bins, G, color=sunset2(0.1), marker='o', markersize=3, lw=1)
    ax.plot(bins, bins, color='black', lw=2, label="Perfect equality")

    ax.legend(fontsize=15, loc='lower right', frameon=False)
    ax.text(0.75, 0.15, r"$\Delta = $" + f"{Gini:.3f}", fontsize=20, ha='center', va='center', transform=ax.transAxes)

    ax.set_ylabel('Occupation', fontsize=16, labelpad=10)

    ax.tick_params(axis='both', which='major', labelsize = 14)

    ax.set_xlabel('Cumulative proportion of grid cells', fontsize=16, labelpad=10)

    # ax.set_title(f"Gini coefficient = {Gini:.3f}", fontsize=16)

    plt.tight_layout()
    plt.show()

def run_3D_tune_positions(parameters3D, show_prints=True):
    if parameters3D["agg"] == 0:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    else:
        filename = f"3D_initial_configurations/3D_neurons_params_{parameters3D['BC_type']}_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_d{parameters3D['d_mean']:.2f}.txt"
    filepath = Path(filename)

    # =========================
    # COMPROBAR SI EXISTE
    # =========================

    if filepath.exists():
        if show_prints:
            print("Configuration already exists:", filename)

    else:
        if show_prints:
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

        if show_prints:
            print(result.stdout)
            print(result.stderr)
        
    return filename

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

def create_network3D(parameters3D, show_prints = True):

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
    if filepath.exists():
        if show_prints:
            print("Configuration already exists:", output_filename)
    else:
        if show_prints:
            print("Creating configuration:", output_filename)

        result = subprocess.run(["./3D_create_network.exe", filename],
        capture_output=True,
        text=True)
        
        if show_prints:
            print(result.stdout)
            print(result.stderr)

    return output_filename

def load_neuron_params_3D(filename, nrows_skiped=15):
    data = np.loadtxt(filename, skiprows=nrows_skiped)
    return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

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
    n_grid = 10

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


def load_A(path):
    A = np.loadtxt(path, dtype=np.uint8)
    np.fill_diagonal(A, 0)
    return A

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

def load_neuron_params_3D(filename, nrows_skiped=15):
    data = np.loadtxt(filename, skiprows=nrows_skiped)
    return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

def load_A(path):
    A = np.loadtxt(path, dtype=np.uint8)
    np.fill_diagonal(A, 0)
    return A

def neuron_map_adj3D(filename_in, filename_out, parameters3D):

    sunset2 = load_cmap('Sunset2', cmap_type='continuous')

    # =========================
    # CARGA
    # =========================
    X, Y, Z, Soma_Diameter, Dendrite_Diameter, Axon_Length = load_neuron_params_3D(filename_in)
    A = load_A(filename_out)

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

def compute_GNA_from_raster(firings_t, firings_i, N, SIM_TIME, W=25, indices_are_1_based=True):
    """
    Calcula la actividad global de red (GNA) como:
        GNA(t) = nº de neuronas activas en la ventana / N

    donde una neurona cuenta solo UNA vez por ventana,
    aunque dispare varios spikes dentro de ella.

    Parámetros
    ----------
    firings_t : array
        Tiempos de disparo (en ms, enteros)
    firings_i : array
        Índices de neurona asociados a cada spike
    N : int
        Número total de neuronas
    SIM_TIME : int
        Tiempo total de simulación (ms)
    W : int
        Tamaño de la ventana deslizante (ms)
    indices_are_1_based : bool
        True si firings_i va de 1..N
        False si firings_i va de 0..N-1

    Devuelve
    --------
    GNA : array de longitud SIM_TIME
        Fracción de neuronas activas en cada instante
    GNA_percent : array de longitud SIM_TIME
        Igual que GNA pero en %
    active_neurons_count : array de longitud SIM_TIME
        Número de neuronas activas en cada instante
    """

    firings_t = np.asarray(firings_t, dtype=int)
    firings_i = np.asarray(firings_i, dtype=int)

    # Pasar índices a 0-based si están en 1..N
    if indices_are_1_based:
        neuron_ids = firings_i - 1
    else:
        neuron_ids = firings_i.copy()

    # Filtrar por seguridad spikes fuera de rango
    valid = (
        (firings_t >= 1) & (firings_t <= SIM_TIME) &
        (neuron_ids >= 0) & (neuron_ids < N)
    )
    firings_t = firings_t[valid]
    neuron_ids = neuron_ids[valid]

    # -----------------------------------------------------
    # 1) Matriz binaria de actividad: activity[n, t] = 1 si
    #    la neurona n dispara en el instante t
    # -----------------------------------------------------
    activity = np.zeros((N, SIM_TIME), dtype=np.uint8)

    # firings_t va de 1..SIM_TIME, así que restamos 1 para indexar 0..SIM_TIME-1
    activity[neuron_ids, firings_t - 1] = 1

    # -----------------------------------------------------
    # 2) Para cada neurona, calculamos si ha estado activa
    #    al menos una vez dentro de la ventana de tamaño W
    #
    #    Usamos convolución temporal por neurona:
    #    si la suma en la ventana > 0, esa neurona está activa
    # -----------------------------------------------------
    kernel = np.ones(W, dtype=int)

    active_in_window = np.zeros_like(activity, dtype=np.uint8)

    for n in range(N):
        counts_in_window = np.convolve(activity[n], kernel, mode='same')
        active_in_window[n] = (counts_in_window > 0).astype(np.uint8)

    # -----------------------------------------------------
    # 3) Número de neuronas activas por instante
    # -----------------------------------------------------
    active_neurons_count = np.sum(active_in_window, axis=0)

    # Fracción y porcentaje
    GNA = active_neurons_count / N
    # GNA_percent = 100 * GNA

    return GNA, active_neurons_count


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


def modify_A_alpha(A, alpha, parameters2D, seed=None):
    """
    A: matriz de conteos (número de intersecciones m_ij)
    alpha: probabilidad de éxito por intento
    parameters2D: diccionario con los parámetros de la red 2D (necesarios para calcular la probabilidad efectiva)
    seed: semilla para reproducibilidad (opcional)
    """

    if not (0 <= alpha <= 1):
        raise ValueError("alpha debe estar entre 0 y 1")

    rng = np.random.default_rng(seed)

    # Probabilidad efectiva por par (i,j)
    P = 1.0 - (1.0 - alpha) ** A

    # Sorteo Bernoulli con esa probabilidad
    A_new = (rng.random(A.shape) < P).astype(np.uint8)

    # Quitar autoconexiones
    np.fill_diagonal(A_new, 0)

    if parameters2D["agg"] == 0:
        output_filename = f"2D_dynamics/networks_for_dynamics/2D_adjacency_matrix_{parameters2D['BC_type']}_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}.txt"
    else:
        output_filename = f"2D_dynamics/networks_for_dynamics/2D_adjacency_matrix_{parameters2D['BC_type']}_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_d{parameters2D['d_mean']:.2f}_a{alpha:.2f}.txt"

    pd.DataFrame(A_new).to_csv(output_filename, sep=' ', index=False, header=False)

    return A_new, output_filename

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