import os

outname = 'myFilelistReCAF.h'
out = open(outname,'w')

out.write("#ifndef %s\n"%(outname.replace('.','_')))
out.write("#define %s\n"%(outname.replace('.','_')))

out.write('static std::vector<std::string> inputFilesReCAF = {\n')
#os.system('ls -1 /pnfs/icarus/scratch/users/jskim/mc/SBNworkshopApril2021__corsika_numu_BNB/CAF/sbncode__v09_27_00_02_FixRangeProton/48388331_*/gen*.root &> tmp_filelist_pnfs.txt')
os.system('ls -1 /pnfs/icarus/scratch/users/jskim/mc/SBNworkshopApril2021__corsika_numu_BNB/FlatCAF/sbncode__v09_27_00_02_FixRangeProton/*.root &> tmp_filelist_pnfs.txt')
os.system('pnfsToXRootD <tmp_filelist_pnfs.txt> tmp_filelist.txt')

lines = open('tmp_filelist.txt').readlines()

for line in lines:
  line = line.strip('\n').replace('root://','"root://').replace('.root','.root",')
  out.write(line+'\n')

out.write('};\n')
out.write('#endif\n')
