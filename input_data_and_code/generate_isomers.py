import os

with open('actives_final.ism') as f: # read in actives_final.ism (as available from DUD-E)
	mol = 1
	for line in f.readlines():
		line = line.split() # for every active, generate a net molecular formula with Open Babel and save the number of the molecule (1-based counting) and its formula in a txt file
		os.system('string=`echo \"' + line[0] + '\" | obabel -i smi -o txt --append formula`; echo ' + str(mol) + ' $string  >> formulae.txt')
		mol += 1
		
num_lines_actives = sum(1 for _ in open('actives_final.ism'))

num_lines_formulae = sum(1 for _ in open('formulae.txt'))

if num_lines_actives != num_lines_formulae:
	print("WARNING!!! DIFFERENT NUMBER OF LINES IN actives_final.ism and formulae.txt IN", os.getcwd(), "!!!")
	
with open('formulae.txt') as f: # for every active from actives_final.ism, read its net formula and do a 10-minute MAYGEN multiprocessor (-m) run to generate isomers. Adjust this, if necessary
	for line in f.readlines():
		line = line.split()
		if not os.path.isfile('../' + line[1] + '.smi'):
			os.system('timeout 10m java -jar /home/user/MAYGEN/target/MAYGEN-1.8.jar -f ' + line[1] + ' -m -v -t -smi -o ../') # adjust the path to your MAYGEN installation
		

