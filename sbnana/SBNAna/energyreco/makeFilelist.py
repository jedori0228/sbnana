import os

wcs = [
'0/*/',
'1/*/',
'2/*/',
'3/*/',
'4/*/',
'5/*/',
'6/*/',
'7/*/',
'8/*/',
#'*/*/',
]

whichDir = 'corsika_nue_BNB'
#whichDir = 'corsika_nueintrinsic_BNB'
#whichDir = 'intimecosmic'

outname = 'myFilelist.h'
outname = 'myFilelistShort.h'

out = open(outname,'w')
out.write('std::vector<std::string> inputFiles = {\n')

for wc in wcs:

  os.system('ls -1 /pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/'+whichDir+'/'+wc+'*.root &> tmp_filelist.txt')

  lines = open('tmp_filelist.txt').readlines()

  for line in lines:
    line = line.strip('\n').replace('/pnfs','"/pnfs').replace('.root','.root",')
    out.write(line+'\n')

out.write('};\n')


