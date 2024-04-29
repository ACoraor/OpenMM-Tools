import os
from subtool import *
from joblib import Parallel, delayed
from clean_thermo import clean_all


def main():
	fns = [ele for ele in os.listdir('.') if "log.lammps." in ele and 
		".all" in ele]
	n = len(fns)
	#Load each log file, process, move on
	curr_dir = os.getcwd()
	#Get basename
	curr_dir = curr_dir[curr_dir[:-1].rfind('/')+1:]
	print("Creating copies at directory %s" % curr_dir)
	lgf = fns
	Parallel(n_jobs=-1)(delayed(make_copy)(fn) for fn in lgf)
	
	#Clean all remaining thermos
	Parallel(n_jobs=-1)(delayed(clean_all)(nr) for nr in range(n))




def make_copy(fname):
	"""Copy the logfile from ../files/fname to inb4(fname,"shuffle"). Extract
	all actual thermo information, and write only the timestep, my_etotal, and temp.
	"""
	print("Copying file %s" % fname)
	writeln = False
	with open(fname,'r') as f:
		with open(inb4(fname,".thermo"),'w') as g: 
			#Write heading
			g.write("#Timestep\tPE\tTemp\n")
			for line in f:
				#Identify line first
				if line[:4] == "Step":
					writeln = True
					continue
				elif line[:4] == "Loop" or line[:4].isalpha():
					writeln = False
					continue

				if writeln:
					#Get correct columns, then write to g.
					cols = str(line).split()
					vals = [cols[0],cols[6],cols[-1]]
					g.write("\t".join(vals) + "\n")

	



if __name__ == '__main__':
	main()
