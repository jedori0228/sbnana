import os

Samples = [
#["MC_NUMI_Nu_Cosmics", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root"],
["MC_NUMI_InTimeCosmics", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root"],
#["DATA_Run7568_NUMI", "/pnfs/icarus/persistent/users/jskim/data/run_7568/flatcaf/v09_45_00/220315_CRTT0Fixed/NUMI/flatcaf_*"],
#["DATA_Run7568_OffbeamNUMI", "/pnfs/icarus/persistent/users/jskim/data/run_7568/flatcaf/v09_45_00/220315_CRTT0Fixed/OffbeamNUMI/flatcaf_*"],
]

for Sample in Samples:

  filelistName = Sample[0]
  filePath = Sample[1]

  outname = 'myFilelist_%s.h'%(filelistName)
  out = open(outname,'w')

  out.write("#ifndef %s\n"%(outname.replace('.','_')))
  out.write("#define %s\n"%(outname.replace('.','_')))

  out.write('static std::vector<std::string> inputFiles_%s = {\n'%(filelistName))

  ## flatcaf, xroodt
  os.system('ls -1 %s &> tmp_filelist_pnfs.txt'%(filePath))

  os.system('pnfsToXRootD <tmp_filelist_pnfs.txt> tmp_filelist.txt')

  lines = open('tmp_filelist.txt').readlines()

  print(filelistName, len(lines))

  for line in lines:

    ## xrootd
    line = line.strip('\n').replace('root://','"root://').replace('.root','.root",')

    ## pnfs
    #line = line.strip('\n').replace('/pnfs','"/pnfs').replace('.root','.root",')

    out.write(line+'\n')

  out.write('};\n')
  out.write('#endif\n')

  os.system('rm tmp_filelist.txt')
  os.system('rm tmp_filelist_pnfs.txt')
