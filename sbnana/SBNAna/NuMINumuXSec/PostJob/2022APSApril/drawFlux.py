## This is not a flux, but a rate of neutrino interaction (CC)

import os,ROOT
from array import array
import tdrstyle
import canvas_margin
from mylib import *

tdrstyle.setTDRStyle()
ROOT.TH1.AddDirectory(False)

ROOT.gROOT.SetBatch(True)

## input

## inputs

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/HistoProducer/"

inputDirName_BeamMC = 'NUMI_Nu_Cosmics/v09_45_00/'

inputDirName = "FMScoreLT12_NuScoreGT0p5_CRTNuMIWindow"

commonFname = "output__Nominal__Incl_PerInt__CutFlow.root"

filepullpath = inputBaseDir+inputDirName+'/'+inputDirName_BeamMC+'/'+commonFname

print('@@ Input sample : %s'%(filepullpath))
f_BeamMC = ROOT.TFile(filepullpath)

outdirBase = "./plots/Flux/"+inputDirName_BeamMC+"/"
os.system("mkdir -p "+outdirBase)

Samples = [
["NuMuCC", "#nu_{#mu} CC", ROOT.kRed],
["NuECC", "#nu_{e} CC", ROOT.kBlue],
]

c = ROOT.TCanvas("c", "", 800, 800)
c.cd()

c.SetLeftMargin(0.20)

h_dummy = ROOT.TH1D("h_dummy", "", 50, 0., 5.)
h_dummy.GetXaxis().SetTitle("Neutrion energy (GeV)")
h_dummy.GetXaxis().SetRangeUser(0, 5)
h_dummy.GetYaxis().SetTitle("Expected number of interactions")
h_dummy.GetYaxis().SetTitleOffset(1.5)
h_dummy.Draw("hist")

latex_POT = ROOT.TLatex()
latex_POT.SetNDC()
latex_POT.SetTextSize(0.035)
latex_POT.DrawLatex(0.70, 0.96, "NuMI, 6e20 POT")

lg = ROOT.TLegend(0.6, 0.8, 0.94, 0.94)
lg.SetBorderSize(0)
lg.SetFillStyle(0)

nRebin = 2
y_max = -999
for Sample in Samples:

  sampleName = Sample[0]
  sampleAlias = Sample[1]
  sampleColor = Sample[2]

  dirName = sampleName+"_NoCut"
  h = f_BeamMC.Get(dirName+"/NeutrinoTruthE_"+dirName)
  h.Rebin(nRebin)
  h.SetLineColor(sampleColor)
  h.SetLineWidth(3)
  h.Draw("histsame")

  lg.AddEntry(h, sampleAlias, "l")

  y_max = max(y_max, h.GetMaximum())

h_dummy.GetYaxis().SetRangeUser(0, 1.1*y_max)
h_dummy.GetYaxis().SetMaxDigits(2)

c.cd()

lg.Draw()

c.SaveAs(outdirBase+"/plot.pdf")

c.SetLogy()
h_dummy.GetYaxis().SetRangeUser(1e2, 50*y_max)
c.SaveAs(outdirBase+"/plot_Logy.pdf")

c.Close()
