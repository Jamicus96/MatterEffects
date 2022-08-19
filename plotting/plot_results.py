### RUN SURVIVAL PROB IN A LOOP AND PLOT ##

import numpy as np
import matplotlib.pyplot as plt
import argparse
import json


def argparser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run matter effect neutrino oscillation in different modes.')
    parser.add_argument('--file', '-f', type=str, dest='file',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Matter_effects/Results/json/MatOsc_stats.json',
                        help='Input file with simulation results.')
    args = parser.parse_args()

    return args
        


def main():
    args = argparser()

    # Read in results
    f = open(args.file)
    data = json.load(f)
    f.close()

    sys_err = 0.01

    L_max = []
    #E_max = np.zeros(len(data))
    P_times = []
    P_times_err = []
    for key in data:
        if data[key]['E_max'] == 1000.0:
            L_max.append(data[key]['L_max'])
            P_times.append(data[key]['time_P']['<x>'])
            P_times_err.append(np.sqrt(data[key]['time_P']['std^2'] + sys_err**2))

    # Plotting
    # plt.plot(L, Pglobesvac, color='red', label='Vacuum Case (GLoBES)')
    # plt.plot(L, Pglobes, color='orange', label='With Matter Effects (GLoBES)')
    # plt.plot(L, Pvac, label='Vacuum Case (my calc)', color='green', linestyle='dashed')
    plt.plot(L_max, P_times, label='With Matter Effects (my calc)', color='blue', linestyle='dashed')
    plt.vlines(L_max, P_times - P_times_err, P_times + P_times_err, color='blue')
    plt.xlabel("L_max (km)")
    plt.legend(loc='best')
    # Title = 'Comparing Survival Probability for {}MeV neutrino,\n\
    #     with Matter Effects (constant density {}g/cm^3) to Vacuum Case'.format(E[0], rho[0])
    # plt.title(Title)

    # # Choose y label (assuming the same for all entries)
    # init_flavour = init_flavour_lst[0]
    # final_flavour = final_flavour_lst[0]
    # if (anti[0] == 1):
    #     if (init_flavour == 0):
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\nu_e \rightarrow \nu_e$ survival probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\nu_e \rightarrow \nu_\mu$ transition probability')
    #         else:
    #             plt.ylabel(r'$\nu_e \rightarrow \nu_\tau$ transition probability')
    #     elif (init_flavour == 1):
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\nu_\mu \rightarrow \nu_e$ transition probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\nu_\mu \rightarrow \nu_\mu$ survival probability')
    #         else:
    #             plt.ylabel(r'$\nu_\mu \rightarrow \nu_\tau$ transition probability')
    #     else:
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\nu_\tau \rightarrow \nu_e$ transition probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\nu_\tau \rightarrow \nu_\mu$ transition probability')
    #         else:
    #             plt.ylabel(r'$\nu_\tau \rightarrow \nu_\tau$ survival probability')
    # else :
    #     if (init_flavour == 0):
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_e$ survival probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_\mu$ transition probability')
    #         else:
    #             plt.ylabel(r'$\overline{\nu}_e \rightarrow \overline{\nu}_\tau$ transition probability')
    #     elif (init_flavour == 1):
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_e$ transition probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_\mu$ survival probability')
    #         else:
    #             plt.ylabel(r'$\overline{\nu}_\mu \rightarrow \overline{\nu}_\tau$ transition probability')
    #     else:
    #         if (final_flavour == 0):
    #             plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_e$ transition probability')
    #         elif (final_flavour == 1):
    #             plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_\mu$ transition probability')
    #         else:
    #             plt.ylabel(r'$\overline{\nu}_\tau \rightarrow \overline{\nu}_\tau$ survival probability')

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