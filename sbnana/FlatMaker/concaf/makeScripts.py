import os

Samples = [

["mc", "NUMI_Nu_Cosmics", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_Nu_Cosmics/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf_*.root", "220520_DefaultRelease_MCP2022A_IcarusProd"],
["mc", "NUMI_in-time_Cosmics2", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics2/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf_*.root", "220520_DefaultRelease_MCP2022A_IcarusProd"],

["data", "run_7897", "v09_37_02_04", "/pnfs/icarus/persistent/users/jskim/data/run_7897/caf/v09_37_02_04/220520_DefaultRelease/NUMI/out/*_*/*stage1.caf.root", "220520_DefaultRelease", "NUMI" ],

]

os.system('mkdir -p grid_dir')

for Sample in Samples:

  sampleType = Sample[0]
  sampleName = Sample[1]
  release = Sample[2]
  sampleCAFPath = Sample[3]
  jobName = Sample[4]
  dataStream = "" if len(Sample)<6 else Sample[5]

  outDir = "/pnfs/icarus/persistent/users/jskim/%s/%s/flatcaf/%s/%s/"%(sampleType,sampleName,release,jobName)
  if sampleType=="data":
    outDir = "/pnfs/icarus/persistent/users/jskim/%s/%s/flatcaf/%s/%s/%s/"%(sampleType,sampleName,release,jobName,dataStream)

  os.system('rm grid_dir/*')

  os.system('ls -1 %s &> tmp.txt'%(sampleCAFPath))
  os.system('pnfsToXRootD <tmp.txt> tmp_xrootd.txt')
  lines = open('tmp_xrootd.txt').readlines()
  NInputfiles = len(lines)

  print("@@ %s, number of files = %d"%(sampleName,len(lines)))

  NTarget = 20

  flistForEachJob = []
  for i in range(0,NTarget):
    flistForEachJob.append( [] )

  for i_line in range(0,len(lines)):
    line = lines[i_line].strip('\n')
    flistForEachJob[i_line%NTarget].append(line)

  for i_flist in range(0,len(flistForEachJob)):

    flist = flistForEachJob[i_flist]
    out = open('grid_dir/run_%s.sh'%(i_flist),'w')
    out.write('#!/bin/bash\n')
    out.write('outDest=$1\n')
    cmd = 'concat_cafs'
    for i_f in range(0,len(flist)):
      out.write('echo "[run_%s.sh] input %d : %s"\n'%(i_flist, i_f, flist[i_f]))
      cmd += ' '+flist[i_f]
    cmd += ' ${outDest}/concat_caf_%d.root\n'%(i_flist)
    cmd += 'flatten_caf ${outDest}/concat_caf_%d.root ${outDest}/flatcaf_%d.root'%(i_flist,i_flist)

    out.write(cmd)
    out.close()

  os.system('tar cf grid_dir.tar grid_dir/')

  submitCMD = '''jobsub_submit \\
-G icarus \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
-l '+SingularityImage=\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\"' \\
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \\
--tar_file_name "dropbox://$(pwd)/grid_dir.tar" \\
--email-to jae.sung.kim.3426@gmail.com \\
-N 20 \\
"file://$(pwd)/grid_executable.sh" \\
"%s"'''%(outDir)

  print(submitCMD)

  os.system(submitCMD)
