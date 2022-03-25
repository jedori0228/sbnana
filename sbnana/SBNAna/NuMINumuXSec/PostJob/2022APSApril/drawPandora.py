import os,ROOT
from array import array
import tdrstyle
import canvas_margin
import mylib
import InTimeCosmicSpillCalculator

ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()
ROOT.TH1.AddDirectory(False)

ROOT.gStyle.SetOptFit(1111)

## Normalization in this plot
POTtoNormalize = 6e20
str_POTtoNormalize = "6e20 POT"

print("@@ POTtoNormalize = %1.2e"%(POTtoNormalize))

## input

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/HistoProducer/"

inputDirName_BeamMC = 'NUMI_Nu_Cosmics/v09_45_00/'

inputDirName = "FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

commonFname = "output__Nominal__Incl_PerInt__CutFlow.root"

f_BeamMC = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_BeamMC+'/'+commonFname)

## Which sample 

Sample = "NuMuCCSignalDef"

cutName = "NoCut"

## Prep

tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextSize(0.045)

## output

outDirBase = './plots/'+inputDirName+'/PandoraPlots/'+inputDirName_BeamMC+'/'
os.system('mkdir -p '+outDirBase)

## Track score

dirName = Sample+"_NoCut"

for ptl in ["Muon", "Proton"]:

  h_TrackScore = f_BeamMC.Get("%s/Truth%sMatchedTrackScore_%s"%(dirName,ptl,dirName))
  h_ShowerScore = f_BeamMC.Get("%s/Truth%sMatchedShowerScore_%s"%(dirName,ptl,dirName))

  c1 = ROOT.TCanvas('c1', '', 800, 800)
  c1.cd()

  h_Score = h_TrackScore.Clone('h_Score')
  h_Score.GetXaxis().SetTitleOffset(1.4)
  h_Score.GetXaxis().SetTitleSize(0.05)
  h_Score.GetXaxis().SetTitle("Track score of %s-matched PF particle"%(ptl.lower()))
  h_Score.GetYaxis().SetTitleOffset(1.5)
  h_Score.GetYaxis().SetTitle("Fraction")
  h_Score.Add(h_ShowerScore)
  h_Score.Scale(1./h_Score.Integral())
  h_Score.Draw("hist")

  trackFrac = h_TrackScore.Integral()/(h_TrackScore.Integral()+h_ShowerScore.Integral())
  tl.DrawLatex(0.2, 0.8, 'Score>0.5 = %1.2f%%'%(trackFrac))

  tl.DrawLatex(0.38, 0.96, "ICARUS Simulation Preliminary")

  c1.SaveAs(outDirBase+"/%sTrackScore.pdf"%(ptl))
  c1.Close()

