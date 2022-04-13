### RUN SURVIVAL PROB IN A LOOP AND PLOT ##

import numpy as np
import argparse
import os


# Optional input arguments
# def argparser():
#     parser = argparse.ArgumentParser(
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description='Run matter effect neutrino oscillation in different modes.')
#     parser.add_argument('--antinu', '-a', type=bool, dest='antinu', choices=[True, False],
#                         default=False, help='True to simulate antineutrinos, False for neutrinos.')
#     args = parser.parse_args()

#     return args


def plot_results(L, E, rho):
    ### IF PLOTTING ###
    import matplotlib.pyplot as plt
    # Read in results
    f = open('results.txt', 'r')
    data = np.loadtxt(f, skiprows=0)

    # neutrino info saved to one number (anti * <init_flavour><final_flavour>, as in <first digit><second digit>)
    anti = np.sign(data[:, 0])
    final_flavour = np.abs(data[:, 0]) % 10
    init_flavour = (np.abs(data[:, 0]) - final_flavour) / 10
    E = data[:, 1]
    rho = data[:, 2]
    L = data[:, 3]
    P = data[:, 0]
    Pvac = data[:, 1]
    Pglobes = data[:, 2]
    Pglobesvac = data[:, 3]

    # Check if all input oscillation is the same
    # result = ((np.max(anti) == np.min(anti)) and (np.max(init_flavour) == np.min(init_flavour)) and (np.max(final_flavour) == np.min(final_flavour)))

    flavours = np.array(['e', '\mu', '\tau'])
    if (anti == 1):
        nu = '\nu_'
    else :
        nu = '\overline{\nu}_'
    y_label = '$' + nu + flavours[int(init_flavour)] + ' \rightarrow ' + nu + flavours[int(final_flavour)] + '$ transition probability'

    plt.plot(L, Pvac, label='Vacuum Case (my calc)')
    plt.plot(L, P, label='With Matter Effects (my calc)', linestyle='dashed')
    plt.plot(L, Pglobesvac, label='Vacuum Case (GLoBES)', linestyle='dashed')
    plt.plot(L, Pglobes, label='With Matter Effects (GLoBES)', linestyle='dashed')
    plt.xlabel("Baseline (km)")
    plt.ylabel(r'$\nu_' + flavours[int(init_flavour)] + ' \rightarrow \nu_' + flavours[int(final_flavour)] + '$ transition probability')
    plt.legend(loc='best')
    Title = 'Comparing Survival Probability for {}MeV neutrino,\n\
        with Matter Effects (constant density {}g/cm^3) to Vacuum Case'.format(E[0], rho[0])
    plt.title(Title)
    plt.show()


def main():

    # Run simulation
    # Compute_oscillations(L, E, rho, args.antinu)

    # Plot results, if desired
    if args.plot:
        plot_results(L, E, rho)


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