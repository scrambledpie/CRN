import subprocess as sp
import numpy as np
import argparse
from datetime import datetime
import os
from fork0_to_csc import CSC_NAMES, callbash, get_hostname
import shutil

home = sp.check_output(['echo $HOME'], shell=True).decode()[:-1]
cpu = get_hostname()


# Call this script from terminal on CSC machine with arguments
#   exp_script: abs path to the script on CSC storage
#   exp_num: int, number of experiments that will be run
#            where exp_script will be called with one integer argument
#   basedir: root folder within which to make new dirname and pass to exp_script
#
#   This script will create a new unique folder, then divide the work up amongst
#   all of the desktops.

def get_git_branch():
    branches = sp.check_output(['git', 'branch']).decode().split("\n")
    for b in branches:
        if "*" in b:
            return b

def host_print(*args):
    """ performs standard print but puts the hostname infront """
    print(cpu + ": ", *args)

def create_exp_dir(exp_script, base_dir="NONAME"):
    """
    Makes a new unique folder inside base_dir, copies source code into the new folder.
    It is assumed that exp_script is within a git repo. The new folder name is appended
    to a txt file full of folder names and numbered by its line number.
    """

    # read the ID.txt file to get the integer uid by counting experiments
    with open( home + '/forkinghellpython/python_savefiles') as f:
        i=0
        for i, _ in enumerate(f):
            pass
    uid = str(i+1)

    # get timestamp
    timestamp = str(datetime.now()).replace(' ', '.')[:-7]

    # git pull, get the commit hash and root dir TODO: strip out \n safe ?
    os.chdir(os.path.dirname(exp_script))
    host_print("current working dir: ", os.getcwd())
    host_print("git branch: " + get_git_branch())
    git_hash = sp.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode()[:-1]
    git_root = sp.check_output(['git', 'rev-parse', '--show-toplevel']).decode()[:-1]

    # read current user
    user = sp.check_output(['whoami']).decode()[:-1]

    # set default base_dir = ~/RESULTS/(repo name)/
    if base_dir == "NONAME":
        repo_name = git_root.split("/")[-1]
        base_dir = home + "/RESULTS/" + repo_name + "/"

        #just make sure the RESULTS dir exists
        if not os.path.exists(home + "/RESULTS"):
            os.makedirs(home + "/RESULTS")

    # construct the new full directory name!
    dirname = base_dir + uid +"."+ user + "."+ timestamp + "."+ git_hash

    # make the dir if it does not exist
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        os.makedirs(dirname+'/StdOut')
        os.makedirs(dirname+'/StdErr')
        os.makedirs(dirname+'/res')
        host_print("Made new dir:  ", dirname)
        host_print("Made new dir:  ", dirname+'/StdOut')
        host_print("Made new dir:  ", dirname+'/StdErr')
        host_print("Made new dir:  ", dirname+'/res')
        host_print("Made new dir:  ", dirname+'/src')


    # append name to ID.txt and create directory
    with open( home + '/forkinghellpython/python_savefiles', 'a') as f:
        f.write(dirname + '\n')
    host_print("Appended name of new root dir to "+ home + "/forkinghellpython/python_savefiles")
    
    # copy source code
    host_print('Backup git repo to '+ dirname + '/src')
    shutil.copytree(git_root, dirname + '/src')
    host_print('Backup complete.')    
    
    return dirname


######################### THE MAIN EVENT ################################

def run(args):

    print("\n")
    host_print("##############  FORK 1 #################")

    N_MACHINES = len(CSC_NAMES)

    # make a new unique dir within the base dir
    dirname = create_exp_dir(args.exp_script, args.basedir)

    # list of jobs and shuffle them
    exp_ids = np.arange(args.exp_num)
    np.random.shuffle(exp_ids)

    # get current conda environment
    conda_env = os.environ['CONDA_DEFAULT_ENV']
    host_print("Conda env: ",conda_env)

    # add verbose flag
    if args.v:
        vflag = " -v"
        host_print("Verbose mode, all StdOut to terminal")
    else:
        vflag = ""
        host_print("Non-verbose mode, all StdOut to", dirname+"/StdOut/")

    # Divide the jobs over CSC machines, give double to kumeta + keiko
    split_ratios = np.ones(len(CSC_NAMES), dtype=float)
    split_ratios += np.array([csc is "kumeta" or csc is "keiko" for csc in CSC_NAMES]) * 0.6
    split_ratios = np.cumsum(split_ratios)
    split_points = np.round(args.exp_num * split_ratios[:-1] / split_ratios[-1]).astype('int')

    # split_points = np.round(np.arange(1, N_MACHINES) * args.exp_num / (N_MACHINES)).astype('int')
    splits = np.split(exp_ids, split_points)
    host_print("Splitting "+str(args.exp_num)+" jobs over ",str(N_MACHINES)," workers.\n\n")

    # For loop over CSC machines
    for name, split in zip(CSC_NAMES, splits):
        # the list of experiment numbers for this CSC machine
        fork_args = [f'--k {i} ' for i in split]

        # call each machine and allocate it some jobs
        # command: python fork_again_on_csc.py (/cscstorage/script.py) (new dirname) (sub list of jobs)
        cmd = "ssh " + name + " 'source $HOME/.bashrc; conda activate "+ conda_env+"; python ~/forkinghellpython/fork2_on_csc.py " + \
                 args.exp_script + " " + dirname + " " + " ".join(fork_args) + vflag + "'&"
        # print(cmd)
        callbash(cmd)

    
    # DONE!


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run experiments on CSC desktops')
    parser.add_argument('exp_script', type=str, help='Experiment script')
    parser.add_argument('exp_num', type=int, help='Number of experiments to run')
    parser.add_argument('--basedir',default="NONAME", type=str, help='Folder to store result in')
    parser.add_argument('-v', action='store_true', help="verbose mode, all output to terminal")
    

    args = parser.parse_args()

    # Just make sure the intended script exists first,
    # let's catch spelling mistakes early!
    if not os.path.isfile(args.exp_script):
        host_print()
        host_print("###############  ERROR  ####################")
        host_print(args.exp_script + " is not a file!!!!! Spelling mistake?")
        import sys; sys.exit()


    # print(args.basedir)
    run(args)
