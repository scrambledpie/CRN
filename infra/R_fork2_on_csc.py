import argparse
import multiprocessing as mp
import os
import socket
import subprocess as sp

from fork0_to_csc import N_PROCESSES, callbash, get_hostname

cpu = get_hostname()


def host_print(*args):
    """ performs standard print but puts the hostname infront """
    print(cpu + ": ", *args)


def run_xp(ARGS):  # k, command_prefix, dirname):
    """ calls the bash command and stdout and stderr are piped into files to reqad later."""
    k = ARGS[0]
    command_prefix = ARGS[1]
    dirname = ARGS[2]
    k = int(k)
    f_out = open(dirname + "/StdOut/" + str(k) + ".txt", 'w')
    f_err = open(dirname + "/StdErr/" + str(k) + ".txt", 'w')
    _ = sp.call(command_prefix + str(k), shell=True, stdout=f_out, stderr=f_err)
    f_out.close()
    f_err.close()
    return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run experiments on CSC desktops')
    parser.add_argument('exp_script', type=str, help='Experiment script')
    parser.add_argument('dirname', type=str, help='Experiment directory')
    parser.add_argument('--k', type=int, action='append', help='Experiment id')
    parser.add_argument('-v', action='store_true', help="verbose mode, all output to terminal")

    args = parser.parse_args()

    # get current conda environment
    conda_env = os.environ['CONDA_DEFAULT_ENV']

    # print out job numbers (tab alignement differs for keiko+kumeta)
    if "keiko" in cpu or "kumeta" in cpu:
        host_print(f'\t\t###### FORK 2 ###### received jobs: {args.k}')
    else:
        host_print(f'\t###### FORK 2 ###### received jobs: {args.k}')

    # parallelize those badboys
    if cpu is "keiko" or cpu is "kumeta":
        pool = mp.Pool(processes=12)
    else:
        pool = mp.Pool(processes=8)

    script_path = os.path.join(args.dirname, 'src') 
    # script_file = os.path.join('scripts', args.exp_script.split('/')[-1])
    script_file = args.exp_script.split('/')[-1]
    command_prefix = "source $HOME/.bashrc; cd " + script_path + "; nice -n 10 Rscript " + script_file + " " + args.dirname + " "

    if args.k is not None:
        # run all the jobs!
        if args.v:
            # verbose mode, output to terminal
            ARGS = [command_prefix + str(k) for k in args.k]
            pool.map(callbash, ARGS)
        else:
            # store StdOut, StdErr in files
            ARGS = [(k, command_prefix, args.dirname) for k in args.k]
            pool.map(run_xp, ARGS)
            # [pool.apply(run_xp, args=(k, command_prefix, args.dirname)) for k in args.k]
