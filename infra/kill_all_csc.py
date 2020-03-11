import subprocess as sp
if __name__ == '__main__':

    def callbash(cmd):
        sp.call(["/bin/bash", "-c", cmd])

    CSC_NAMES = ["rilyeno", "torta", "adobo", "bulalo", "kinilaw", "okoy",
                 "embutido", "jamon", "caldereta", "dinuguan", "lechon",
                 "niliga", "inihaw", "halabos", "sinuglaw", "keiko", "kumeta"]


    for name in CSC_NAMES:
        callbash("ssh " + name + " 'killall python'&")