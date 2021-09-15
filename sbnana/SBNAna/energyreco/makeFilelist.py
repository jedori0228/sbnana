import os

wcs = [
'*/*/',
]

whichDir = 'corsika_nue_BNB'
#whichDir = 'corsika_nueintrinsic_BNB'
#whichDir = 'intimecosmic'

outname = 'myFilelist.h'
#outname = 'myFilelistShort.h'

out = open(outname,'w')
out.write("#ifndef %s\n"%(outname.replace('.','_')))
out.write("#define %s\n"%(outname.replace('.','_')))
out.write('static std::vector<std::string> inputFiles = {\n')

for wc in wcs:

  os.system('ls -1 /pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/'+whichDir+'/'+wc+'*.root &> tmp_filelist_pnfs.txt')
  os.system('pnfsToXRootD <tmp_filelist_pnfs.txt> tmp_filelist.txt')

  lines = open('tmp_filelist.txt').readlines()

  for line in lines:
    line = line.strip('\n').replace('root://','"root://').replace('.root','.root",')
    out.write(line+'\n')

out.write('};\n')
out.write('#endif\n')
