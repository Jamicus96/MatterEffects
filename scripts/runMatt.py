from math import degrees
import numpy as np
import argparse
import subprocess
import os
import time
import json

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--example_job', '-ej', type=str, dest='example_job',
                        default='job_scripts/jobArray.job', help='Example job script.')
    parser.add_argument('--info_repo', '-ir', type=str, dest='info_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Matter_effects/Results/job_info/', help='Folder to save job info text files in.')
    parser.add_argument('--stats_repo', '-sr', type=str, dest='stats_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Matter_effects/Results/output_stats/', help='Folder to save timing result text files in.')
    parser.add_argument('--json_repo', '-jr', type=str, dest='json_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Matter_effects/Results/json/', help='Folder to save final json stats in.')

    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                    default=True, help='print and save extra info')
    parser.add_argument('---step', '-s', type=str, dest='step', default='all', choices=['sim', 'stats', 'all'],
                        help='which step of the process is it in?')

    args = parser.parse_args()
    return args


### MISCELLANEOUS FUNCTIONS ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('scripts/runMatt.py')]
    if repo_address == '':
        firt_char = None
    else:
        firt_char = repo_address[0]
    if firt_char != '/':
        repo_address = os.getcwd() + '/' + repo_address
    return repo_address

def checkRepo(repo_address, verbose=False):
    '''Check format of repo address is right to be used here. Also if repo does not exist, create it.'''

    # Format address
    new_address = repo_address
    if new_address == '.':
        new_address = ''
    elif new_address != '':
        if (new_address[int(len(new_address) - 1)] != '/'):
            new_address += '/'
        
        # If directory does not exists, make it
        if not os.path.isdir(new_address):
            os.mkdir(new_address)
            if verbose:
                print('Created new directory: ', new_address)

    return new_address

def makeJobArrayScript(example_jobScript, save_job_folder, commandList_address, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = save_job_folder + 'MatOsc.job'

    output_logFile_address = save_job_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_MatOsc.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif '/Address/CommandList.txt' in line:
            new_line = line.replace('/Address/CommandList.txt', commandList_address, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def checkJobsDone(jobName_str, wait_time, verbose):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    # Turns out the name of the job is only the 10 first characters of the job file name

    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE, shell=True).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
                if len(jobName_str) > 10:
                    jobName_str = jobName_str[:9]
                if jobName_str in line:
                    running = True
                    if verbose:
                        print('Waiting for jobs to finish...')
                    break
        time.sleep(wait_time)

    return True

### SIMS FUNCTIONS ###

def setSimVals(save_stats_folder):
    # Number of times to repeat the same measurement
    repeat_measurement = 10

    # Set unchanging values
    init_flavour = '0'    # 0 = e, 1 = mu, 2 = tau
    final_flavour = '0'   # 0 = e, 1 = mu, 2 = tau
    anti = '-1'           # -1 =  antineutrino, 1 =  neutrino
    print_probs = '0'     # 1  =  true, 0  =  false

    L_N = '100'           # number of data points between L_min and L_max
    E_N = '100'           # number of data points between E_min and E_max
    rho_N = '100'         # number of data points between rho_min and rho_max
    L_log = '0'           # 1  =  true, 0  =  false
    E_log = '0'           # 1  =  true, 0  =  false
    rho_log = '0'         # 1  =  true, 0  =  false
    L_min = '0.01'        # km
    E_min = '0.5'         # MeV
    rho_min = '0'         # g/cm^3

    # Set changing values
    L_max_list  = np.linspace(0.1, 1000, 10)    # km
    E_max_list = np.linspace(0.5, 1000, 10)     # MeV
    rho_max_list = np.linspace(2.7, 2.7, 10)     # g/cm^3

    info = []
    for L_max in L_max_list:
        L_max = str(L_max)
        for E_max in E_max_list:
            E_max = str(E_max)
            for rho_max in rho_max_list:
                rho_max = str(rho_max)
                for i in range(repeat_measurement):
                    i = str(i)
                    output_file_address = save_stats_folder + 'results_' +  L_min  + '_' +  L_max  + '_' +  L_N  + '_' +  L_log \
                                            + '_' +  E_min  + '_' +  E_max  + '_' +  E_N  + '_' +  E_log \
                                            + '_' +  rho_min  + '_' +  rho_max  + '_' +  rho_N  + '_' +  rho_log \
                                            + '_' +  init_flavour  + '_' +  final_flavour  + '_' +  anti  + '_' +  print_probs \
                                            + '_' +  i  + '.txt'
                    info_line = [L_min, L_max, L_N, L_log, E_min, E_max, E_N, E_log, rho_min, rho_max, rho_N, rho_log,\
                                init_flavour, final_flavour, anti, output_file_address, print_probs]
                    info.append(info_line)

    return info

def runSims(args):
    '''Runs simulations based in input information'''
    print('Running runSims().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    with open(repo_address + args.example_job, "r") as f:
        example_jobScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)
    save_job_folder = checkRepo(args.info_repo, args.verbose)
    
    ### MAKE JOB SCRIPT TO RUN THE SIMULATIONS ###
    print('Creatingand job script and info file...')

    # Creat info file, tha contains all the commands
    input_info = setSimVals(save_stats_folder)
    commandList_address = save_job_folder + 'command_list.txt'
    commandList_file = open(commandList_address, 'w')
    command_base = repo_address + 'scripts/Matter_Oscillations.exe'
    for info_line in input_info:
        command = command_base
        for info in info_line:
            command += ' ' + str(info)
        commandList_file.write(command + '\n')
    commandList_file.close()

    # Create job file that will execute all these commands
    job_address = makeJobArrayScript(example_jobScript, save_job_folder, commandList_address, args.verbose)


    ### RUN JOB SCRIPTS ###
    print('Submitting job array...')
    command = 'qsub -t 1-' + str(len(input_info)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True 

### Analysis functions ###

def rolling_av(old_av, new_x, n):
    '''n is the number of values used to compute this new average (n = 1, 2, 3, ...)'''
    if n == 1:
        return new_x
    else:
        return old_av + ((new_x - old_av) / n)

def rolling_var(old_var, new_x, new_av, n):
    '''n is the number of values used to compute this new variance (n = 1, 2, 3, ...)'''
    if n == 1:
        return 0.0
    else:
        return ((n - 1.0) / n) * old_var + ((1.0 / n) + (1 / (n - 1.0)**2)) * (new_x - new_av)**2

def listDir(dir):
    '''List all the files in a directory'''
    init_list = os.listdir(dir)
    file_list = []
    for item in init_list:
        if os.path.isfile(dir + item):
            file_list.append(item)
    return file_list

def check_entry_exists(info, temp_info, name_info_format):
    line = 0
    for info_line in info:
        exists = True
        for stats_info_name in name_info_format:
            if info_line[stats_info_name] != temp_info[stats_info_name]:
                exists = False
                break
        if exists:
            return True, line
        line += 1
    return False, -1
                
def getInfo(save_stats_folder, files):

    info = []
    name_info_format = ['L_min', 'L_max', 'L_N', 'L_log', 'E_min', 'E_max', 'E_N', 'E_log',\
                        'rho_min', 'rho_max', 'rho_N', 'rho_log', 'init_flavour', 'final_flavour',\
                        'anti', 'print_probs']
    stats_info_format = ['time_P', 'time_P_vac', 'time_P_globes', 'time_P_vac_globes', 'time_P_specific']
    for file in files:
        # Get info in name
        file_info = file.replace('.txt', '').split('_')
        if file_info[0] != 'results':
            continue
        file_info = file_info[1:]

        # Create new entry
        temp_info = {}
        good_entry = False
        for i in range(len(name_info_format)):
            temp_info[name_info_format[i]] = float(file_info[i])

        # Get file info
        file = open(save_stats_folder + file, 'r')
        lines = file.readlines()
        file.close()
        for line in lines:
            if len(line) > 12:
                if line[:12] == '- - - - - - ':
                    stats_info = line.replace('- - - - - - ', '').split(' ')
                    for i in range(len(stats_info)):
                        temp_info[stats_info_format[i]] = {}
                        temp_info[stats_info_format[i]]['<x>'] = float(stats_info[i])
                        if i == len(stats_info_format) - 1:
                            good_entry = True
                    break

        if good_entry:
            # Check if this entry already exists
            exists, line_num = check_entry_exists(info, temp_info, name_info_format)
            
            if not exists:  # If it doesn't, just add to info
                for stats_info_name in stats_info_format: # Set std^2 = 0 (only one val so far)
                    temp_info[stats_info_name]['std^2'] = 0.0
                temp_info['n_data'] = 1.0 # First data point
                # Add to stats
                info.append(temp_info)

            else: # If it does, update values with rolling average and std^2
                info[line_num]['n_data'] += 1.0 # New data point we're adding
                for stats_info_name in stats_info_format: # Set std^2 = 0 (only one val so far)
                    info[line_num][stats_info_name]['<x>'] = rolling_av(info[line_num][stats_info_name]['<x>'], temp_info[stats_info_name]['<x>'], info[line_num]['n_data'])
                    info[line_num][stats_info_name]['std^2'] = np.abs(rolling_var(info[line_num][stats_info_name]['std^2'], temp_info[stats_info_name]['<x>'],
                                                                        info[line_num][stats_info_name]['<x>'], info[line_num]['n_data']))

    return info

def readStats(args):
    # Make sure folders are of the correct format to  use later
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)
    save_json_folder = checkRepo(args.json_repo, args.verbose)

    # Wait for jobs to finish
    checkJobsDone('MatOsc.job', 10, args.verbose)

    # List files in repo
    files = listDir(save_stats_folder)
    print(files)

    # Get all relevent info from file name
    info = getInfo(save_stats_folder, files)
    print(info)

    # Write info to json file
    table = {}
    i = 0
    for line in info:
        table[str(i)] = line
        i += 1
    with open(save_json_folder + 'MatOsc_stats.json', 'w') as f:
        json.dump(table, f)


### OTHER ###

def runAll(args):
    runSims(args)
    readStats(args)
    return True

### MAIN ###

def main():
    # read in argument
    args = argparser()

    work_modes = {
        'sim': runSims,
        'stats': readStats,
        'all': runAll
    }

    result = work_modes[args.step](args)

if __name__ == '__main__':
    main()