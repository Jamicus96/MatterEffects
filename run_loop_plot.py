### RUN SURVIVAL PROB IN A LOOP AND PLOT ##

import numpy as np
import argparse
import os


# Oscillation constants
theta12 = np.arcsin(np.sqrt(0.307))
theta13 = np.arcsin(np.sqrt(2.18e-2))
theta23 = np.arcsin(np.sqrt(0.545))     # Normal Hierarchy
# theta23 = np.arcsin(np.sqrt(0.547))   # Inverted Hierarchy
delta = 1.36 * np.pi                    # Not sure about the hierarchy
m21 = 7.53e-5                           # (Delta m_21^2, in eV^2)
m31 = 0.0025283                         # (Delta m_31^2, in eV^2)
# m32 = m31 - m21                       # (Delta m_32^2, in eV^2)


# Optional input arguments
def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run matter effect neutrino oscillation in different modes.')
    parser.add_argument('--antinu', '-a', type=bool, dest='antinu', choices=[True, False],
                        default=False, help='True to simulate antineutrinos, False for neutrinos.')
    parser.add_argument('--plot', '-p', type=bool, dest='plot', choices=[True, False],
                        default=False, help="Plot results directly? Doesn't work on hpc")
    args = parser.parse_args()

    return args


def print_to_file(filename, info):
    # If file exists, delete and recreate it
    file_exists = os.path.exists(filename)
    if file_exists:
        os.system('rm ' + filename)
    f = open(filename, 'w')

    # Write to file line by line
    if isinstance(info, list) or isinstance(info, np.ndarray):
        for i in range(info):
            f.write(info[i])
    else:
        f.write(info)


def calculate_PMNS():
    # Compute useful constants
    s12 = np.sin(theta12)
    s13 = np.sin(theta13)
    s23 = np.sin(theta23)
    c12 = np.cos(theta12)
    c13 = np.cos(theta13)
    c23 = np.cos(theta23)
    c_delta = np.cos(delta)
    s_delta = np.sin(delta)

    # Real part of PMNS matrix
    U_r = np.zeros((3, 3))
    U_r[0, 0] = c12 * c13;                                  U_r[0, 1] = s12 * c13;                                  U_r[0, 2] = s13 * c_delta
    U_r[1, 0] = - s12 * c23 - s13 * s23 * c12 * c_delta;    U_r[1, 1] = c12 * c23 - s12 * s13 * s23 * c_delta;      U_r[1, 2] = s23 * c13
    U_r[2, 0] = s12 * s23 - s13 * c12 * c23 * c_delta;      U_r[2, 1] = - s23 * c12 - s12 * s13 * c23 * c_delta;    U_r[2, 2] = c13 * c23

    # Imaginary part of PMNS matrix
    U_i = np.zeros((3, 3))
    U_i[0, 0] = 0.0;                                        U_i[0, 1] = 0.0;                                        U_i[0, 2] = - s13 * s_delta
    U_i[1, 0] = - s13 * s23 * c12 * s_delta;                U_i[1, 1] = - s12 * s13 * s23 * s_delta;                U_i[1, 2] = 0.0
    U_i[2, 0] = - s13 * c12 * c23 * s_delta;                U_i[2, 1] = - s12 * s13 * c23 * s_delta;                U_i[2, 2] = 0.0

    return U_r, U_i


def Compute_oscillations(L, E, rho, anti):
    # Create new file (delete first, if it exists)
    file_exists = os.path.exists('results.txt')
    if file_exists:
        os.system('rm results.txt')
    open('results.txt', 'w')

    ### Run simulation for each parameter ###
    for i in range(len(L)):
        os.system('./Matter_Oscillations.exe {} {} {} 1'.format(E, L[i], rho))


def plot_results(L, E, rho):
    ### IF PLOTTING ###
    import matplotlib.pyplot as plt
    # Read in results
    f = open('results.txt', 'r')
    data = np.loadtxt(f,skiprows=0)

    P = data[:,0]
    Pvac = data[:,1]
    Pglobes = data[:,2]
    Pvacglobes = data[:,3]

    plt.plot(L, Pvac, label='Vacuum Case (my calc)')
    plt.plot(L, P, label='With Matter Effects (my calc)', linestyle='dashed')
    plt.plot(L, Pvacglobes, label='Vacuum Case (GLoBES)', linestyle='dashed')
    plt.plot(L, Pglobes, label='With Matter Effects (GLoBES)', linestyle='dashed')
    plt.xlabel("Baseline (km)")
    plt.ylabel(r"$\nu_\mu \rightarrow \nu_e$ transition probability")
    plt.legend(loc='best')
    Title = 'Comparing Survival Probability for {}MeV neutrino,\n\
        with Matter Effects (constant density {}g/cm^3) to Vacuum Case'.format(E, rho)
    plt.title(Title)
    plt.show()


def main():
    # Get arguments
    args = argparser()

    # Print oscillation constants to file, to be used by C++ code, ensuring consistency
    osc_consts = [theta12, theta13, theta23, delta, m21, m31]
    print_to_file('oscillation_constants.txt', osc_consts)

    # Compute PMNS matrix elements and print them to file (for vacuum oscillation calculation)
    U_r, U_i = calculate_PMNS()
    PMNS_elements = []
    for i in range(3):
        for j in range(3):
            PMNS_elements.append(U_r[i][j])
    for i in range(3):
        for j in range(3):
            PMNS_elements.append(U_i[i][j])
    print_to_file('PMNS_m2_data.txt', PMNS_elements)

    # Compute constants needed for algorithm, outside of loop (delete file these are printed to first)
    filename = 'Constants.txt'
    file_exists = os.path.exists(filename)
    if file_exists:
        os.system('rm ' + filename)
    os.system('./Mat_Os_Consts.exe')

    ### Simulation parameters ###
    N = 500
    L = np.linspace(0.0, 100.0, N)
    E = 1.0
    rho_lithosphere = 2.7  #g/cm^3
    rho = 0.0  #g/cm^3

    # Run simulation
    Compute_oscillations(L, E, rho, args.antinu)

    # Plot results, if desired
    if args.plot:
        plot_results(L, E, rho)





#  Convert matter density to electron density using 
#  converstion factors from GLOBES-3.0.11/src/glb_probability.h:
#  GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
#  GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#  Also note:
#  V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below which vacuum algorithms are used */
#  Anyway, used to convert from matter density rho (g/cm^3) to matter potential V (glb_probability.c),
#  Where: V = sqrt(2)*G_F*N_e, and is given in eV

# Ne = rho / 4.821e-15 #electron density cm^-3
# Ne *= 1e6 #(cm^-3 to m^3)