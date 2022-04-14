### RUN SURVIVAL PROB IN A LOOP AND PLOT ##

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run matter effect neutrino oscillation in different modes.')
    parser.add_argument('--file', '-f', type=bool, dest='file',
                        default='results.txt', help='Input file with simulation results.')
    args = parser.parse_args()

    return args
        


def main():
    args = argparser()

    # Read in results
    f = open(args.file, 'r')
    data = np.loadtxt(f, skiprows=0)

    # neutrino info saved to one number (<anti><init_flavour><final_flavour>, as in <first digit><second digit><third digit>)
    anti = np.sign(data[:, 0]).astype(int)
    final_flavour_lst = np.abs(data[:, 0]).astype(int) % 10
    init_flavour_lst = ((np.abs(data[:, 0]) - final_flavour_lst) / 10).astype(int) % 10
    E = data[:, 1]
    rho = data[:, 2]
    L = data[:, 3]
    P = data[:, 4]
    Pvac = data[:, 5]
    Pglobes = data[:, 6]
    Pglobesvac = data[:, 7]

    # Check if all input oscillation is the same
    # result = ((np.max(anti) == np.min(anti)) and (np.max(init_flavour) == np.min(init_flavour)) and (np.max(final_flavour) == np.min(final_flavour)))

    # Plotting
    plt.plot(L, Pvac, label='Vacuum Case (my calc)')
    plt.plot(L, P, label='With Matter Effects (my calc)', linestyle='dashed')
    plt.plot(L, Pglobesvac, label='Vacuum Case (GLoBES)', linestyle='dashed')
    plt.plot(L, Pglobes, label='With Matter Effects (GLoBES)', linestyle='dashed')
    plt.xlabel("Baseline (km)")
    plt.legend(loc='best')
    Title = 'Comparing Survival Probability for {}MeV neutrino,\n\
        with Matter Effects (constant density {}g/cm^3) to Vacuum Case'.format(E[0], rho[0])
    plt.title(Title)

    # Choose y label (assuming the same for all entries)
    init_flavour = init_flavour_lst[0]
    final_flavour = final_flavour_lst[0]
    print(anti[0])
    print(init_flavour)
    print(final_flavour)
    if (anti[0] == 1):
        if (init_flavour == 0):
            if (final_flavour == 0):
                plt.ylabel(r'$\nu_e \rightarrow \nu_e$ survival probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\nu_e \rightarrow \nu_\mu$ transition probability')
            else:
                plt.ylabel(r'$\nu_e \rightarrow \nu_\tau$ transition probability')
        elif (init_flavour == 1):
            if (final_flavour == 0):
                plt.ylabel(r'$\nu_\mu \rightarrow \nu_e$ transition probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\nu_\mu \rightarrow \nu_\mu$ survival probability')
            else:
                plt.ylabel(r'$\nu_\mu \rightarrow \nu_\tau$ transition probability')
        else:
            if (final_flavour == 0):
                plt.ylabel(r'$\nu_\tau \rightarrow \nu_e$ transition probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\nu_\tau \rightarrow \nu_\mu$ transition probability')
            else:
                plt.ylabel(r'$\nu_\tau \rightarrow \nu_\tau$ survival probability')
    else :
        if (init_flavour == 0):
            if (final_flavour == 0):
                plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_e$ survival probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_\mu$ transition probability')
            else:
                plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_\tau$ transition probability')
        elif (init_flavour == 1):
            if (final_flavour == 0):
                plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_e$ transition probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_\mu$ survival probability')
            else:
                plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_\tau$ transition probability')
        else:
            if (final_flavour == 0):
                plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_e$ transition probability')
            elif (final_flavour == 1):
                plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_\mu$ transition probability')
            else:
                plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_\tau$ survival probability')

    plt.show()

if __name__ == '__main__':
    main()


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