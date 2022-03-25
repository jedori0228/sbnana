import os,ROOT
from array import array
import tdrstyle
import canvas_margin
from mylib import *

ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()
ROOT.TH1.AddDirectory(False)

ROOT.gStyle.SetOptFit(1111)

TargetPOT = 6e20
str_TargetPOT = "6e20 POT"
latex_POT = ROOT.TLatex()
latex_POT.SetNDC()
latex_POT.SetTextSize(0.035)

## inputs

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/QuickEfficiency/"

inputCats = dict()
#inputCats["Default"] = "FMScoreRemoved_NuScoreRemoved_DefaultCRTHitVeto"
#inputCats["NUMIWindow"] = "FMScoreRemoved_NuScoreRemoved_NUMIWindowCRTHitVeto"
#inputCats["NUMIWindowPE150"] = "FMScoreRemoved_NuScoreRemoved_NUMIWindowPE150CRTHitVeto"

inputCats["Default"] = "FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

commonFname = "output.root"

inputDirNameInfos = [

#["v09_28_01_01_01", "NUMI_Nu_Cosmics", "Default"],
#["v09_28_01_01_01", "NUMI_in-time_Cosmics", "Default"],

#["v09_37_01_03p01", "NUMI_Nu_Cosmics", "Default"],
#["v09_37_01_03p01", "NUMI_in-time_Cosmics", "Default"],

["v09_45_00", "NUMI_Nu_Cosmics", "Default"],
["v09_45_00", "NUMI_in-time_Cosmics", "Default"],

#["v09_45_00", "NUMI_Nu_Cosmics", "NUMIWindow"],
#["v09_45_00", "NUMI_in-time_Cosmics", "NUMIWindow"],

#["v09_45_00", "NUMI_Nu_Cosmics", "NUMIWindowPE150"],
#["v09_45_00", "NUMI_in-time_Cosmics", "NUMIWindowPE150"],

]

TFs = dict()
for inputDirNameInfo in inputDirNameInfos:
  key = inputDirNameInfo[0]+"/"+inputDirNameInfo[1]+"/"+inputDirNameInfo[2]
  TFs[key] = ROOT.TFile(inputBaseDir+inputCats[inputDirNameInfo[2]]+'/'+inputDirNameInfo[1]+'/'+inputDirNameInfo[0]+'/'+commonFname)

## samples

samples = [

#[TFs['v09_28_01_01_01/NUMI_Nu_Cosmics/Default'], 'NuMuCC', '#nu#mu CC (Old)'],
#[TFs['v09_28_01_01_01/NUMI_Nu_Cosmics/Default'], 'Cosmic', 'Out-time cosmic (Old)'],
#[TFs['v09_28_01_01_01/NUMI_in-time_Cosmics/Default'], 'AllSamples', 'In-time cosmic (Old)'],

#[TFs['v09_37_01_03p01/NUMI_Nu_Cosmics/Default'], 'NuMuCC', '#nu#mu CC (Prod.)'],
#[TFs['v09_37_01_03p01/NUMI_Nu_Cosmics/Default'], 'Cosmic', 'Out-time cosmic (Prod.)'],
#[TFs['v09_37_01_03p01/NUMI_in-time_Cosmics/Default'], 'AllSamples', 'In-time cosmic (Prod.)'],

[TFs['v09_45_00/NUMI_Nu_Cosmics/Default'], 'NuMuCC', '#nu#mu CC (v0945)'],
#[TFs['v09_45_00/NUMI_Nu_Cosmics/Default'], 'Cosmic', 'Out-time cosmic (v0945)'],
[TFs['v09_45_00/NUMI_in-time_Cosmics/Default'], 'AllSamples', 'In-time cosmic (v0945)'],

#[TFs['v09_45_00/NUMI_Nu_Cosmics/NUMIWindow'], 'NuMuCC', '#nu#mu CC (v0945, Test)'],
#[TFs['v09_45_00/NUMI_Nu_Cosmics/NUMIWindow'], 'Cosmic', 'Out-time cosmic (v0945, Test)'],
#[TFs['v09_45_00/NUMI_in-time_Cosmics/NUMIWindow'], 'AllSamples', 'In-time cosmic (v0945, Test)'],

#[TFs['v09_45_00/NUMI_Nu_Cosmics/NUMIWindowPE150'], 'NuMuCC', '#nu#mu CC (v0945, Test2)'],
#[TFs['v09_45_00/NUMI_Nu_Cosmics/NUMIWindowPE150'], 'Cosmic', 'Out-time cosmic (v0945, Test2)'],
#[TFs['v09_45_00/NUMI_in-time_Cosmics/NUMIWindowPE150'], 'AllSamples', 'In-time cosmic (v0945, Test2)'],

]


cutFlows = [

#['BaseSelection', ['NoCut', 'RecoFiducial', 'NotClearCosmic']],
['CRTFromNoCut', ['NoCut', 'NoCut_SideCRTVeto', 'NoCut_CRTVeto']],
#['CRTFromRecoFiducial', ['RecoFiducial', 'RecoFiducial_SideCRTVeto', 'RecoFiducial_CRTVeto']],
#['CRTFromNotClearCosmic', ['NotClearCosmic', 'NotClearCosmic_SideCRTVeto', 'NotClearCosmic_CRTVeto']],
#['CRTFromFMScoreCut', ['FMScoreCut', 'FMScoreCut_SideCRTVeto', 'FMScoreCut_CRTVeto']],

]

varName = 'CountSlice'

for cutFlow in cutFlows:

  setName = cutFlow[0]
  cuts = cutFlow[1]

  Header = setName
  for sample in samples:
    sampleName = sample[2]
    Header += '\t'+sampleName
  print(Header)

  y_Start = []

  for i_cut in range(0,len(cuts)):

    cut = cuts[i_cut]

    out = cut

    for i_sample in range(0,len(samples)):

      sample = samples[i_sample]

      sampleName = sample[1]

      h = sample[0].Get(sampleName+"_"+cut+"/CountSlice_"+sampleName+"_"+cut)
      y = h.GetBinContent(1)

      eff = 1.
      if i_cut==0:
        y_Start.append(y)
      else:
        eff = y/y_Start[i_sample]

      out += '\t%1.3f'%(eff)

    print(out)
