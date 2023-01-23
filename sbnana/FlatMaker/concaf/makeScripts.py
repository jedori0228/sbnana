import os

Samples = [

#["mc", "NUMI_Nu_Cosmics", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_Nu_Cosmics/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf_*.root", "220520_DefaultRelease_MCP2022A_IcarusProd"],
#["mc", "NUMI_in-time_Cosmics2", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics2/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf_*.root", "220520_DefaultRelease_MCP2022A_IcarusProd"],

#["mc", "BNB_Nu_Cosmics", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/BNB_Nu_Cosmics/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf_*.root", "220520_DefaultRelease_MCP2022A_IcarusProd"],
#["mc", "BNB_Nu_Cosmics", "v09_37_02_04", "/pnfs/icarus/persistent/users/jskim/mc/BNB_Nu_Cosmics/caf/v09_37_02_04/220524_ProdNuSyst/out/33020419_*/prod*.root", "220525_NuSyst_MCP2022A_IcarusProd"],
#["mc", "BNB_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/BNB_Nu_Cosmics_Ovb/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_Default/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_Default"],
#["mc", "BNB_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/BNB_Nu_Cosmics_Ovb/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_PanCautiousCRMode/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_PanCautiousCRMode"],

#["mc", "NUMI_Nu_Cosmics", "v09_37_02_04", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_Nu_Cosmics/mc/reconstructed/icaruscode/v09_37_02_04/caf/*/*/caf*.root", "icarus_numi_nu_cosmics_v09_37_02_04_caf"],
#["mc", "NUMI_Nu_Cosmics", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_Default/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_Default"],
#["mc", "NUMI_Nu_Cosmics", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_PanCautiousCRMode_MaxNuCos9999/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_PanCautiousCRMode_MaxNuCos9999"],
#["mc", "NUMI_Nu_Cosmics", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_PanCautiousCRMode/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_PanCautiousCRMode"],

#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/220913_Default/out/38678230_*/prod*.root", "icarus_numi_nu_cosmics_ovb_v09_37_02_04_caf"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/220727_NuSyst_ExtendedRR/out/36446710_*/prodcorsika*.root", "220727_NuSyst_ExtendedRR"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_Default/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_Default"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/221017_PandoraCosmicRecoTest_PanCautiousCRMode/out/*_*/prod*.root", "221017_PandoraCosmicRecoTest_PanCautiousCRMode"],

#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221027_PandoraCosmicVertexTest_Default/out/*_*/prod*.root", "TEST_221027_PandoraCosmicVertexTest_Default"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221027_PandoraCosmicVertexTest_DaughterSet-9999/out/*_*/prod*.root", "TEST_221027_PandoraCosmicVertexTest_DaughterSet-9999"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221027_PandoraCosmicVertexTest_ParentSet-9999/out/*_*/prod*.root", "TEST_221027_PandoraCosmicVertexTest_ParentSet-9999"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221101_PandoraCosmicVertexTest_CAFNotSkipSlice/out/*_*/prod*.root", "TEST_221101_PandoraCosmicVertexTest_CAFNotSkipSlice"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221102_PandoraCosmicVertexTest_VtxDefaultY-20/out/*_*/prod*.root", "TEST_221102_PandoraCosmicVertexTest_VtxDefaultY-20"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221102_PandoraCosmicVertexTest_BugFix/out/*_*/prod*.root", "TEST_221102_PandoraCosmicVertexTest_BugFix"],
#["mc", "NUMI_Nu_Cosmics_Ovb", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics_Ovb/caf/v09_37_02_09/TEST_221102_PandoraCosmicVertexTest_BugFix/out/*_*/prod*.root", "TEST_221102_PandoraCosmicVertexTest_BugFix"],

#["mc", "NUMI_in-time_Cosmics_withOverburden2", "v09_37_02_07", "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics_withOverburden2/mc/reconstructed/icaruscode/v09_37_02_07/caf/*/*/caf_*.root", "IcarusProd_2022A_NUMI_in-time_Cosmics_withOverburden_v09_37_02_07_caf"],
#["mc", "NUMI_in-time_Cosmics_withOverburden2", "v09_61_00_01", "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics_withOverburden2/caf/v09_61_00_01/221114_UpdatedMCCalibConst/out/40572999_*/prod*.root", "221114_UpdatedMCCalibConst"],

#["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/220727_ExtendedRR/NUMIMAJORITY/out/58561773_*/data*.root", "220727_ExtendedRR", "NUMIMAJORITY" ],
#["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/220909_NuSyst_ExtendedRR_FMTimeFixedForData/NUMIMAJORITY/out/60777842_*/data*.root", "220909_NuSyst_ExtendedRR_FMTimeFixedForData", "NUMIMAJORITY" ],
#["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/220915_NuSyst_ExtendedRR_FMTimeFixedForData4usShift/NUMIMAJORITY/out/61056822_*/data*.root", "220915_NuSyst_ExtendedRR_FMTimeFixedForData4usShift", "NUMIMAJORITY" ],
["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/221120_PMTCRTMatchingEff/NUMIMAJORITY/out/*_*/data*.root", "221120_PMTCRTMatchingEff", "NUMIMAJORITY" ],
["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/221120_PMTCRTMatchingEff/OffBeamNUMIMAJORITY/out/*_*/data*.root", "221120_PMTCRTMatchingEff", "OffBeamNUMIMAJORITY" ],
["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/221120_PMTCRTMatchingEff/BNBMAJORITY/out/*_*/data*.root", "221120_PMTCRTMatchingEff", "BNBMAJORITY" ],
["data", "run_8515", "v09_37_02_09", "/pnfs/icarus/persistent/users/jskim/data/run_8515/caf/v09_37_02_09/221120_PMTCRTMatchingEff/OffBeamBNBMAJORITY/out/*_*/data*.root", "221120_PMTCRTMatchingEff", "OffBeamBNBMAJORITY" ],

#["data", "run_8437", "v09_37_02_09", "/pnfs/icarus/scratch/users/jlarkin/cafs/numi/8437/data_*.OKTOLOOK.caf.root", "220927_BlindCAF", "NUMIMAJORITY"],
#["data", "run_8439", "v09_37_02_09", "/pnfs/icarus/scratch/users/jlarkin/cafs/numi/8439/data_*.OKTOLOOK.caf.root", "220927_BlindCAF", "NUMIMAJORITY"],
#["data", "run_8446", "v09_37_02_09", "/pnfs/icarus/scratch/users/jlarkin/cafs/numi/8446/data_*.OKTOLOOK.caf.root", "220927_BlindCAF", "NUMIMAJORITY"],

#["data", "run_Run1", "v09_60_01", "/pnfs/icarus/scratch/users/gputnam/DP2022P/run1-majorityD/62375300_*/data_*.OKTOLOOK.caf.root", "221104_Gray_D", "NUMIMAJORITY"],

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
-N %d \\
"file://$(pwd)/grid_executable.sh" \\
"%s"'''%(NTarget,outDir)

  print(submitCMD)

  os.system(submitCMD)
