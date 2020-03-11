import argparse
import os
import subprocess as sp
from time import sleep


def get_hostname():
    return sp.check_output(['hostname'], shell=True).decode()[:-1]


cpu = get_hostname()


def host_print(*args):
    """ performs standard print but puts the hostname infront """
    print(cpu + ": ", *args)


def get_home():
    return sp.check_output(['echo $HOME'], shell=True).decode()[:-1]


N_PROCESSES = 8

#########################################################################################    
#
#          HOW TO USE THIS SCRIPT!!!
#
#########################################################################################
#
# REQUIREMENTS
#   - ssh keys for passwordless login between local->CSC and CSC<->CSC
#   - conda to be installed on CSC account eg ":~$ source $HOME/.bashrc; conda activate TFenv; python myTFscript.py"
#   - tmux installed (sudo apt install tmux) for htop viewer
#
#########################################################################################
#
# TLDR:
# STEP 1. check machines are awake + free and delete bad ones from list below (literally edit the list in this file).
# From local/CSC terminal type "python fork0_to_csc.py; tmux attach -t Nhtop"
# with each call to this fork0, this file is copied to CSC so the list of machines in this file
# is the list of machines that will be used for jobs.
#
# STEP 2: to run the example script 10 times, type the following command in local terminal:
# "python fork0_to_csc.py \$HOME/cond_bayes_opt/scripts/demo_infra_usage.py 10 -v"
# which will call
#  "python (CSC_HOME_DIR)/cond_bayes_opt/scripts/demo_infra_usage.py (new bespoke dirname) (k)"
# for k in 0,...,9 spread over all the CSC machines in the list below. See demo_infra_usage.py
# for further usage. It is also possible to set a custom conda env, and custom output dir, branch.
# the -v flag is for verbose printing stdout+stderr to terminal otherwise it saves to files.
# there is a --nopull flag mean CSC repo will not pull from github.
# 
#########################################################################################
#
# CALLING THIS SCRIPT WITHOUT ARGUMENTS 
#
# just loads a tmux session with loads of htops, 
# use this first to check which CSC machines are free. In bash, type 
#           "python fork0_to_csc.py; tmux attach -t Nhtop".
# and if some machines do not respond or are already in use, delete them from the CSC_NAMES list below.
# This script also copies this file (with the latest CSC_NAMES) over to CSC storage so the next
# fork1_on_csc.py and fork2_on_csc.py executed from CSC storage will have the same copy of CSC_NAMES.
#
#########################################################################################
#
# CALLING THIS SCRIPT WITH ARGUMENTS
#
#   exp_script: absolute path to the script on CSC storage. This script must be callable
#               with the given conda env as "python (exp_script) (dirname) (k)"
#               Therefore treat exp_script as a single master runner with all experiment settings
#               in a lookup table. Then calling this script with arguments (dirname) (k) 
#               executes the experiment with row k of settings and saves output in (dirname).
# 
#   exp_num: int, the total number of experiments that will be run
#            with exp_script with k = 0,1,2,3,...,exp_num
#
#   --basedir: output base dir on csc storage within which a new unique output dirname will be created.
#              repo contents are coptied into dirname and dirname is then passed to 
#              exp_script ( and can be used for saving outputs).
#              Default basedir is $HOME/RESULTS/(repo name)/
#
#   --first_fork: str, (optional) default=jamon, ssh name of first CSC desktop to connect to
#
#   --conda: str, (optional) default="base". the name of the conda environment to be used on the CSC
#
#   --branch: str, (optional) name of git repo branch to use in the repo on CSC 
#             note that (branch setting is persistent so no need to set every time)
#
#   -v: verbose mode, stdout+stderr go to terminal. Otherwise stdout+stderr will be saved to unique .txt files.
#
#   -nopull: this flag STOPS the default behaviour, which is to logon to CSC and "cd (path_exp_script); git pull".
#
#   This fork0_to_csc will copy this fork0_to_csc.py file over to csc so that CSC_MACHINES list 
#   on the CSC storage is
#   up-to-date. fork0_to_csc will then logon to args.first_fork, and call fork1_on_csc. The fork1_on_csc will
#   navigate to basedir, make a new unique dirname, then partition the 1,...,exp_num jobs into sets for each CSC machine
#   in the list above (used for tmux Nhtop). Then fork1_on_csc.py will logon to each csc machine
#   and call the next fork2_on_csc. This next script will then use the provided conda environment and the partition of jobs
#   to execute "python exp_script (dirname) (k)"" for all k in this csc machine partition of 1,..,exp_num jobs.
#   The new dirname is appended to a list of dirnames in godzilla:~/forkinghellpython/python_savefiles.
#
#   The new dirname is saved into godzilla:~/forkinghellpython/python_savefiles
#
# e.g. call from local machine
#    python fork0_to_csc.py /home/maths/phrnaj/MCBO/run_MCBO_exp.py 1000 --first_fork adobo --basedir /home/maths/phrnaj/MCBO_results/ --conda TFgpu
#
# This will result in runs "python home/maths/phrnaj/MCBO/run_MCBO_exp.py basedir  k" repeatedly for k in 1,...,exp_num spread over CSC desktops.
#
##########################################################################################


############################### DEFINE COMPUTER ARRAY #######################################

# set the computers you want to use here, tmux will load to show if they are active.
ALL_CSC_NAMES = ["rilyeno", "torta", "adobo", "bulalo", "kinilaw", "okoy",
                 "embutido", "jamon", "caldereta", "dinuguan", "lechon",
                 "niliga", "inihaw", "halabos", "sinuglaw", "keiko", "kumeta"]

# default list uses all computers, but some may need to be removed.
# working with names is a bitch, instead use numbers (tmux panes).
U = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
U = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

CSC_NAMES = [ALL_CSC_NAMES[i] for i in U]


################################ UTILITY FUNCTIONS ##########################################

# Make a tmux session full of htop windows
def Nhtop(names=CSC_NAMES):
    sp.Popen(['/bin/bash', '-c', 'tmux new -d -s Nhtop'])
    sp.Popen(['/bin/bash', '-c', 'tmux new-window -t Nhtop'])
    sp.Popen(['/bin/bash', '-c', 'tmux swap-window -t Nhtop:0 -s Nhtop:1'])
    sp.Popen(['/bin/bash', '-c', 'tmux kill-window -t Nhtop:1'])

    sleep(0.1)

    # create panes and load htop
    for i in names:
        if i != names[0]:
            sp.call(["/bin/bash", "-c", "tmux split-window -v -t Nhtop"])
            sp.call(["/bin/bash", "-c", "tmux select-layout -t Nhtop tiled"])

        sp.call(["/bin/bash", "-c", "tmux send-keys -t Nhtop 'ssh " + i + "' 'C-m'"])

        # if there are ssh errors from host to CSC, delete ~/.ssh/known_hosts and set this to True
        if 1==11:
            # if there are ssh errors from CSC to CSC, delete known_hosts in CSC storage and uncomment next 2 lines
            # sleep(2)
            # sp.call(["/bin/bash", "-c", "tmux send-keys -t Nhtop 'ssh " + i + "' 'C-m'"])
            sleep(2)
            sp.call(["/bin/bash", "-c", "tmux send-keys -t Nhtop 'yes' 'C-m'"])
        else:
            sp.call(["/bin/bash", "-c", "tmux send-keys -t Nhtop 'htop' 'C-m'"])

        sleep(0.15)


def callbash(cmd, silent=False):
    if silent:
        _ = sp.check_output(["/bin/bash", "-c", cmd])
    else:
        _ = sp.call(["/bin/bash", "-c", cmd], stderr=sp.STDOUT)


################################ THE MAIN EVENT ##########################################

if __name__ == '__main__':
    print("\n")
    print("###### FORK 0 (LOCAL) ######")
    print("htop array loading, type: tmux attach -t Nhtop")
    Nhtop()

    # copy over these forkings files to CSC storage so that 'CSC_NAMES' list is up-to-date.
    print("Copying local fork*.py files to godzilla:~/forkinghellpython/")
    try:
        # make a folder to store these files on CSC
        callbash("ssh godzilla 'mkdir ~/forkinghellpython'")
    except:
        pass
    callbash("scp fork*.py godzilla:~/forkinghellpython/")
    callbash("ssh godzilla 'touch forkinghellpython/python_savefiles'")
    print("Copied")

    # command line arguments
    parser = argparse.ArgumentParser(description='Run experiments on CSC desktops')
    parser.add_argument('exp_script', type=str, help='absolute path on CSC storage to experiment script')
    parser.add_argument('exp_num', type=int, help='Number of experiments to run')
    parser.add_argument('--basedir', type=str, default='NONAME',
                        help='Dest. folder on CSC to backup source code, also passed to exp_script')
    parser.add_argument('--first_fork', type=str, default='jamon', help='which CSC desktop to start from')
    parser.add_argument('--conda', type=str, default='base', help='the conda environment on CSC for exp_script')
    parser.add_argument('--branch', type=str, default=None, help='the git branch to use')
    parser.add_argument('-v', action='store_true', help="verbose mode, all output to terminal")
    parser.add_argument('-nopull', action='store_true', help="Do not perform git pull on CSC")

    # This will crash if there are no arguments
    args = parser.parse_args()
    git_dir = os.path.dirname(args.exp_script)

    # make sure the exp_script is a git repo (fork1 will break otherwise)
    try:
        cmd = "ssh godzilla 'cd " + git_dir + "; git rev-parse --git-dir'"
        AA = sp.check_output(cmd, shell=True)
    except:
        print("The given script file is not in a git repo!")
        import sys;

        sys.exit()

    # login to godzilla and pull repo unless user says no
    if not args.nopull:
        print("Doing git pull on CSC: " + git_dir)
        callbash("ssh godzilla 'cd " + git_dir + "; git pull'")

    # log into godzilla and set the branch, it is persistent after logout/login
    if args.branch is not None:
        cmd = "ssh godzilla 'cd " + git_dir + "; git checkout " + args.branch + "'"
        callbash(cmd)
        print("On CSC, set branch to " + args.branch)

    # define verbose flag
    if args.v:
        vflag = " -v"
    else:
        vflag = ""

    # This will not run if there are no arguments
    print("Calling all CSC machines! Your country needs you!   Over to ..... " + args.first_fork, "\n\n")
    callbash(
        f"ssh {args.first_fork} 'source ~/.bashrc; conda activate {args.conda}; cd ~/forkinghellpython/; python fork1_on_csc.py " + \
        args.exp_script + " " + str(args.exp_num) + " --basedir " + args.basedir + vflag + "'&")
