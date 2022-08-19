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
                        default='/Users/jp643/Documents/Studies/PhD/Antinu/Matter_Effects/Results/MatOsc_stats.json',
                        help='Input file with simulation results.')
    args = parser.parse_args()

    return args
        
def check_match(set_vals, data):
    for name, val in set_vals:
        if val != data[name]:
            return False
    return True


def main():
    args = argparser()

    # Read in results
    f = open(args.file)
    data = json.load(f)
    f.close()

    sys_err = 0.01

    # Set values
    set_vals = [("L_min", 0.01),
                # ("L_N", 100.0),
                # ("L_log", 0.0),
                ("E_min", 0.5),
                # ("E_N", 100.0),
                # ("E_log", 0.0),
                ("rho_min", 0.0),
                # ("rho_N", 100.0),
                # ("rho_log", 0.0),
                ("init_flavour", 0.0),
                ("final_flavour", 0.0),
                ("anti", -1.0),
                
                ("L_max", 1000.0),
                # ("E_max", 1000.0),
                ("rho_max", 2.7)
                ]

    # Ranging values
    running_vals = [
                    # ('L_max', [], 'L_max (km)'),
                    ('E_max', [], 'E_max (MeV)'),
                    # ('rho_max', [], 'rho_max (g/cm^3)')
                    ]

    CPU_times = [('time_P', [], [], 'New algorithm, general flavour', 'blue'),
                ('time_P_specific', [], [], 'New algorithm, specific flavour', 'green'),
                ('time_P_vac', [], [], 'Vacuum calculation using PMNS matrix', 'red'),
                ('time_P_globes', [], [], 'GLoBES algorithm', 'yellow'),
                ('time_P_vac_globes', [], [], 'GLoBES vacuum algorithm', 'orange')]  # (time name, time, time err, plot label, plot colour)

    # Get data wanted
    for key in data:
        if check_match(set_vals, data[key]):
            for x_name, x_vals, x_label in running_vals:
                x_vals.append(data[key][x_name])
            for name, times, time_errs, label, col in CPU_times:
                print(name)
                times.append(data[key][name]['<x>'])
                time_errs.append(np.sqrt(data[key][name]['std^2'] + sys_err**2))

    # Repackage
    np_arr_running_vals = []
    for x_name, x_vals, x_label in running_vals:
        np_arr_running_vals.append((name, np.array(x_vals), x_label))
    np_arr_CPU_times = []
    for name, times, time_errs, label, col in CPU_times:
        np_arr_CPU_times.append((name, np.array(times), np.array(time_errs), label, col))
    # Plotting
    for x_name, x_vals, x_label in np_arr_running_vals:
        for name, times, time_errs, label, col in np_arr_CPU_times:
            plt.scatter(x_vals, times, label=label, color=col)
            plt.vlines(x_vals, times - time_errs, times + time_errs, color=col)

        plt.ylabel('Simulation time (s)')
        plt.xlabel(x_label)
        plt.legend(loc='best')
        plt.title('Comparing Speed of Neutrino Oscillation Algorithms')

        plt.show()

if __name__ == '__main__':
    main()

