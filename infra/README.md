     
  
 # HOW TO USE

 # REQUIREMENTS
  - ssh keys for passwordless login between local->CSC and CSC<->CSC
  - conda to be installed on CSC account eg ":~$ source $HOME/.bashrc; conda activate TFenv; python myTFscript.py"
  - tmux installed locally, (sudo apt install tmux) for htop viewer

########################################################################################

# SEE fork0_to_csc.py for more documentation!

Here is the TLDR:
STEP 1. to check machines that are free and delete from list if necessary,
from local/CSC terminal type "python fork0_to_csc.py; tmux attach -t Nhtop" 

STEP 2: to run a script 1000 times across the csc machines, type the following command:
"python fork0_to_csc.py /cscstorage/myscript.py 1000 /cscstorage/resultsfolder --first_fork adobo --conda myenv"
which will run the script "/cscstorage/myscript.py" stored on CSC storage 1000 times 
using virtual environment myenv, using adobo as the entry point to CSC, and inside /cscstorage/resultsfolder a new UNIQUE
folder is made and git repo src code is copied over.

fork0_to_csc.py will do a git pull on the CSC copy of the repo.
