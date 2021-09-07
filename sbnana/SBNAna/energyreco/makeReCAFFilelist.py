import os

outname = 'myFilelistReCAF.txt'
out = open(outname,'w')
out.write('std::vector<std::string> inputFilesReCAF = {\n')
os.system('ls -1 /pnfs/icarus/scratch/users/jskim/mc/SBNworkshopApril2021__corsika_numu_BNB/CAF/sbncode__v09_27_00_02/47977478_*/gen*.root &> tmp_filelist.txt')

lines = open('tmp_filelist.txt').readlines()

for line in lines:
	line = line.strip('\n').replace('/pnfs','"/pnfs').replace('.root','.root",')
	out.write(line+'\n')

out.write('};\n')

