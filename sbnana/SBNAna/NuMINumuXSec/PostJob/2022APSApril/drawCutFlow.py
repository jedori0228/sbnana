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

latex_POT = ROOT.TLatex()
latex_POT.SetNDC()
latex_POT.SetTextSize(0.035)

## input

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/HistoProducer/"

inputDirName_BeamMC = 'NUMI_Nu_Cosmics/v09_45_00/'
inputDirName_IntimeCosmicMC = 'NUMI_in-time_Cosmics/v09_45_00/'

inputDirName = "FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

commonFname = "output__Nominal__Incl_PerInt__CutFlow.root"

f_BeamMC = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_BeamMC+'/'+commonFname)
f_IntimeCosmicMC = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_IntimeCosmicMC+'/'+commonFname)

nominalIntensity = 48.09116545e12

## For the in-time cosmic MC normalization
m_InTimeCosmicSpillCalculator = InTimeCosmicSpillCalculator.InTimeCosmicSpillCalculator(POTtoNormalize, nominalIntensity)
InTimeCosmicSpill = m_InTimeCosmicSpillCalculator.GetInTimeCosmicSpill(f_BeamMC, f_IntimeCosmicMC)

## Which sample 

Samples = [
[f_BeamMC, 'NuMuCCSignalDef', '#nu_{#mu}-CC', ROOT.kRed, "signal"],
[f_BeamMC, 'NuMuNC', '#nu_{#mu}-NC', ROOT.kBlue, "bkgd"],
[f_BeamMC, 'Cosmic', 'Out-time cosmic', ROOT.kYellow+2, "bkgd"],
[f_IntimeCosmicMC, 'AllSamples', 'In-time cosmic', ROOT.kOrange, "bkgd"],
]

## Efficiency plots first

cutNames = [

['NoCut','No cut'],
['RecoFiducial','Fiducial volume'],
['NotClearCosmic','Not a clear cosmic'],
['FMScore','Flash matching score<12'],
['SideCRTVeto','Side CRT veto'],
['CRTVeto','Top CRT veto'],
['HasMuon','Muon track'],

]

varName = 'CountSlice'

## output

outDirBase = './plots/'+inputDirName+'/Cutflow/'+inputDirName_BeamMC+'/'
os.system('mkdir -p '+outDirBase)

## Draw

c1 = ROOT.TCanvas('c1', '', 1200, 800)
c1.SetBottomMargin(0.20)
c1.cd()

h_dummy = ROOT.TH1D('h_dummy', '', len(cutNames), 0., 1.*len(cutNames))
h_dummy.GetYaxis().SetTitle("Efficiency")
h_dummy.Draw("axis")

dict_h = dict()
for Sample in Samples:
  dict_h[Sample[2]] = ROOT.TH1D('Eff_'+Sample[2], '', len(cutNames), 0., 1.*len(cutNames))
  dict_h[Sample[2]].SetLineColor(Sample[3])

for i in range(0,len(cutNames)):
  h_dummy.GetXaxis().SetBinLabel(i+1, cutNames[i][1])
  for Sample in Samples:
    h = Sample[0].Get('%s_%s/%s_%s_%s'%(Sample[1],cutNames[i][0],varName,Sample[1],cutNames[i][0]))
    dict_h[Sample[2]].SetBinContent(i+1, h.GetBinContent(1))

lg_Eff = ROOT.TLegend(0.25, 0.3, 0.45, 0.50)
lg_Eff.SetBorderSize(0)
lg_Eff.SetFillStyle(0)

grs = []
map_Scaled = dict()
for Sample in Samples:

  ## Eff

  h_Eff = dict_h[Sample[2]].Clone( dict_h[Sample[2]].GetName()+"_Eff" )
  h_Eff.Scale(1./h_Eff.GetBinContent(1))

  gr = mylib.convertToGraph(h_Eff)
  gr.SetMarkerSize(0)
  gr.SetMarkerColor(0)
  gr.SetLineColor(Sample[3])
  gr.SetLineWidth(2)
  gr.Draw("ezsame")

  lg_Eff.AddEntry(gr, Sample[2], 'l')

  grs.append(gr)

  ## Scaledity

  h_Scaled = dict_h[Sample[2]].Clone( dict_h[Sample[2]].GetName()+"_Scaledity" )
  if Sample[2]=='In-time cosmic':

    h_Livetime = Sample[0].Get('%s_%s/%s_%s_%s'%(Sample[1],cutNames[i][0],'Livetime',Sample[1],cutNames[i][0]))
    this_Livetime = h_Livetime.GetBinContent(1)
    print("In-time cosmic MC contains %1.2e spills"%(this_Livetime))
    ## InTime cosmic is already scaled by IntimeCosmicMC_TargetPOT/1e18
    ## https://github.com/SBNSoftware/sbnana/blob/e030b862847d8e84dc35006906276ca58be396ca/sbnana/CAFAna/Core/Spectrum.cxx#L442
    SpillScale = InTimeCosmicSpill/this_Livetime
    print("-> to match InTimeCosmicSpill (%1.2e), In-time MC should be scaled by %1.2e"%(InTimeCosmicSpill, SpillScale))
    alreadyScaled = m_InTimeCosmicSpillCalculator.IntimeCosmicMC_TargetPOT/1e18
    print("-> In-time MC has been already scaled by %1.2e/1e18 = %1.2e"%(m_InTimeCosmicSpillCalculator.IntimeCosmicMC_TargetPOT, alreadyScaled))
    h_Scaled.Scale(SpillScale/alreadyScaled)
    print("In-time cosmic scaled by %1.2e"%(SpillScale/alreadyScaled))
  else:
    h_Scaled.Scale(POTtoNormalize/m_InTimeCosmicSpillCalculator.BeamMC_TargetPOT)
    print("Beam MC scaled by %1.2e"%(POTtoNormalize/m_InTimeCosmicSpillCalculator.BeamMC_TargetPOT))

  map_Scaled[Sample[2]] = h_Scaled.Clone()

hs_Purity = ROOT.THStack("hs_Purity", "")
for i in range(0,len(cutNames)):
  y_SumThisBin = 0.
  for i_Sample in range(0,len(Samples)):
    Sample = Samples[ len(Samples)-i_Sample-1 ]
    y_SumThisBin += map_Scaled[Sample[2]].GetBinContent(i+1)
  for i_Sample in range(0,len(Samples)):
    Sample = Samples[ len(Samples)-i_Sample-1 ]
    map_Scaled[Sample[2]].SetBinContent(i+1, map_Scaled[Sample[2]].GetBinContent(i+1)/y_SumThisBin)

lg_Eff.Draw()

h_dummy.GetYaxis().SetRangeUser(1E-4,2)
c1.SetLogy()

latex_POT.DrawLatex(0.65, 0.96, "ICARUS Simulation Preliminary")

c1.SaveAs(outDirBase+'/Cutflow_Eff.pdf')

c1.cd()
hs_Purity = ROOT.THStack("hs_Purity", "")
for i in range(0,len(cutNames)):
  y_SumThisBin = 0.
  for i_Sample in range(0,len(Samples)):
    Sample = Samples[ len(Samples)-i_Sample-1 ]
    y_SumThisBin += map_Scaled[Sample[2]].GetBinContent(i+1)
  for i_Sample in range(0,len(Samples)):
    Sample = Samples[ len(Samples)-i_Sample-1 ]
    newBinCont = map_Scaled[Sample[2]].GetBinContent(i+1)/y_SumThisBin
    print(i,Sample[2],newBinCont)
    map_Scaled[Sample[2]].SetBinContent(i+1, map_Scaled[Sample[2]].GetBinContent(i+1)/y_SumThisBin)


lg_Pur = ROOT.TLegend(0.25, 0.6, 0.45, 0.80)
lg_Pur.SetBorderSize(0)
lg_Pur.SetFillStyle(1001)
lg_Pur.SetFillColor(ROOT.kWhite)

for i_Sample in range(0,len(Samples)):

  #Sample = Samples[ len(Samples)-i_Sample-1 ]
  Sample = Samples[i_Sample]

  map_Scaled[Sample[2]].SetFillColor(Sample[3])
  hs_Purity.Add(map_Scaled[Sample[2]])
  #lg_Pur.AddEntry(map_Scaled[Sample[2]], Sample[2], "f")

for i_Sample in range(0,len(Samples)):

  Sample = Samples[i_Sample]
  Sample = Samples[ len(Samples)-i_Sample-1 ]

  lg_Pur.AddEntry(map_Scaled[Sample[2]], Sample[2], "f")

h_dummy.Draw("hist")
h_dummy.GetYaxis().SetRangeUser(0., 1.1)
c1.SetLogy(False)
hs_Purity.Draw("histsame")
h_dummy.Draw("axissame")
h_dummy.GetYaxis().SetTitle("Purity")

lg_Pur.Draw()

latex_POT.DrawLatex(0.65, 0.96, "ICARUS Simulation Preliminary")

c1.SaveAs(outDirBase+'/Cutflow_Purity.pdf')

c1.Close()
