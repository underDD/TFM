import subprocess
import numpy as np
from pathlib import Path
from pypalettes import load_cmap
import matplotlib.pyplot as plt

def load_neuron_params(filename, nrows_skiped=13):
        data = np.loadtxt(filename, skiprows=nrows_skiped)
        return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

def run_2D_tune_positions(parameters2D, show_prints=True):

    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}.txt"

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
            str(parameters2D["base_sigma"])
        ], capture_output=True, text=True)

        if show_prints:
            print(result.stdout)
            print(result.stderr)

    return filename

def create_network2D(parameters2D, alphas, show_prints = True):
    
    if parameters2D["agg"] == 0:
        filename = f"2D_initial_configurations/2D_neurons_params_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}.txt"
    else:
        filename = f"2D_initial_configurations/2D_neurons_params_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}.txt"
        
    # --------------------------------

    output_filenames = []
    for a in alphas:
            if parameters2D["agg"] == 0:
                output_filenames.append(f"2D_created_networks/2D_adjacency_matrix_random_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_a{a:.3f}.txt")
            else:
                output_filenames.append(f"2D_created_networks/2D_adjacency_matrix_agg_nc{parameters2D['n_centers']}_s{parameters2D['base_sigma']:.4f}_L{parameters2D['L']:.1f}_rho{parameters2D['rho']:.0f}_l{parameters2D['l_mean']:.2f}_a{a:.3f}.txt")

    output_filepaths = []

    for out_filename in output_filenames:
        filepath = Path(out_filename)
        output_filepaths.append(filepath)

    # print("Output filepaths:", output_filepaths)

    for filepath, a in zip(output_filepaths, alphas):
        if filepath.exists():
            if show_prints:
                print("Configuration already exists:", filepath)
        else:
            if show_prints:
                print("Creating configuration:", filepath)

            result = subprocess.run(["./2D_create_network.exe", f"{a:.3f}", filename],
            capture_output=True,
            text=True)
            
            if show_prints:
                print(output_filenames)

                print(result.stdout)
                print(result.stderr)
        
    return output_filenames

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
        filename = f"3D_initial_configurations/3D_neurons_params_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}.txt"
    else:
        filename = f"3D_initial_configurations/3D_neurons_params_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}.txt"
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
            str(parameters3D["base_sigma"])
        ], capture_output=True, text=True)

        if show_prints:
            print(result.stdout)
            print(result.stderr)
        
    return filename

def create_network3D(parameters3D, alphas, show_prints = True):

    if parameters3D["agg"] == 0:
        filename = f"3D_initial_configurations/3D_neurons_params_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}.txt"
    else:
        filename = f"3D_initial_configurations/3D_neurons_params_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}.txt"

    # --------------------------------

    output_filenames = []
    for a in alphas:
            if parameters3D["agg"] == 0:
                output_filenames.append(f"3D_created_networks/3D_adjacency_matrix_random_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_a{a:.3f}.txt")
            else:
                output_filenames.append(f"3D_created_networks/3D_adjacency_matrix_agg_nc{parameters3D['n_centers']}_s{parameters3D['base_sigma']:.4f}_L{parameters3D['L']:.1f}_rho{parameters3D['rho']:.0f}_l{parameters3D['l_mean']:.2f}_a{a:.3f}.txt")

    output_filepaths = []

    for out_filename in output_filenames:
        filepath = Path(out_filename)
        output_filepaths.append(filepath)

    # print("Output filepaths:", output_filepaths)

    for filepath, a in zip(output_filepaths, alphas):
        if filepath.exists():
            if show_prints:
                print("Configuration already exists:", filepath)
        else:
            if show_prints:
                print("Creating configuration:", filepath)

            result = subprocess.run(["./3D_create_network.exe", f"{a:.3f}", filename],
            capture_output=True,
            text=True)
            if show_prints:
                print(result.stdout)
                print(result.stderr)
    
    return output_filenames

def load_neuron_params_3D(filename, nrows_skiped=14):
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