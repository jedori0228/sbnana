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

lt = ROOT.TLatex()
lt.SetNDC()
lt.SetTextSize(0.045)

## input

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/HistoProducer/"

inputDirName_BeamMC = 'NUMI_Nu_Cosmics/v09_45_00/'

inputCat = "FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

commonFname = "output__Nominal__Incl_PerInt__CutFlow.root"

filepullpath = inputBaseDir+inputCat+'/'+inputDirName_BeamMC+'/'+commonFname

print('@@ Input sample : %s'%(filepullpath))
f_BeamMC = ROOT.TFile(filepullpath)

## Which sample 

sample = "NuMuCCSignalDef"

## Efficiency plots first

varNames = [
'NeutrinoTruthE',
'MuonTruthP',
'MuonTruthNuMICosineTheta',
]

h_empty_legend = ROOT.TH1D("h_empty_legend", "", 1, 0., 1.)

CutSets = [

[
'NoCut',
'RecoFiducial',
],
[
'NoCut',
'RecoFiducial',
'NotClearCosmic',
'FMScore',
],
[
'NoCut',
'RecoFiducial',
'NotClearCosmic',
'FMScore',
'SideCRTVeto',
'CRTVeto',
],
[
'NoCut',
'RecoFiducial',
'NotClearCosmic',
'FMScore',
'SideCRTVeto',
'CRTVeto',
'HasMuon',
],

]

for cutNames in CutSets:

  #outDirBase = './plots/'+inputCat+'/Efficiency/'+inputDirName_BeamMC+'/'
  outDirBase = './plots/'+inputCat+'/Efficiency/'+inputDirName_BeamMC+'/Cut%d/'%(len(cutNames)-1)
  os.system('mkdir -p '+outDirBase)

  for varName in varNames:

    ## Efficiency based on true disbtribution

    histName_Truth = varName
    xTitle = ''
    xMin = 0.
    xMax = 5.
    dx = 0.1
    nRebin = 1
    if varName=='NeutrinoTruthE':
      histName_Truth = 'NeutrinoTruthE'
      xTitle = 'Neutrino energy (GeV)'
      xMax = 10.
      nRebin = 2
    elif varName=='MuonTruthP':
      histName_Truth = 'MuonTruthP'
      xTitle = 'Muon momentum (GeV)'
      nRebin = 2
    elif varName=='MuonTruthCosineTheta':
      xTitle = 'Muon Cos(#theta)'
      xMin = -1.
      xMax = 1.
      dx = 0.1
    elif varName=='MuonTruthNuMICosineTheta':
      xTitle = 'Muon Cos(#theta_{NuMI})'
      xMin = -1.
      xMax = 1.
      dx = 0.1
    elif varName=='ProtonTruthP':
      xTitle = 'Proton momentum (GeV)'

    hs_Truth = []
    colors = [
      ROOT.kBlack,
      ROOT.kRed,
      ROOT.kOrange,
      ROOT.kViolet,
      ROOT.kBlue,
      ROOT.kGreen,
      ROOT.kGray+2,
      ROOT.kGray,
      
    ]
    for i_dir in range(0,len(cutNames)):
      cutName = cutNames[i_dir]
      dirName = sample+'_'+cutName
      tDir = f_BeamMC.Get( dirName )

      h = tDir.Get( histName_Truth+'_'+dirName )
      print(i_dir,cutName)
      hs_Truth.append(h)

    c1 = ROOT.TCanvas('c1', '', 800, 800)

    c1_up = ROOT.TPad("c1_up", "", 0, 0.25, 1, 1)
    c1_down = ROOT.TPad("c1_down", "", 0, 0, 1, 0.25)

    c1, c1_up, c1_down = canvas_margin.canvas_margin(c1, c1_up, c1_down)
    c1.Draw()
    c1_up.Draw()
    c1_down.Draw()

    c1_up.cd()

    h_dummy_up = ROOT.TH1D('h_dummy_up', '', int( (xMax-xMin)/dx ), xMin, xMax)
    h_dummy_up.Draw("axis")

    y_max = -999.

    for i_h in range(0,len(hs_Truth)):
      h = hs_Truth[i_h]
      h.Rebin(nRebin)
      cutName = cutNames[i_h]
      color = colors[i_h]
      h.SetLineColor(color)
      h.Draw("histsame")
      y_max = max(y_max, h.GetMaximum())

    h_dummy_up.GetYaxis().SetRangeUser(1, 1.4*y_max)

    lg_1 = ROOT.TLegend(0.18, 0.69, 0.95, 0.91)
    lg_1.SetNColumns(2)
    lg_1.SetBorderSize(0)
    lg_1.SetFillStyle(0)
    lg_1.AddEntry(hs_Truth[0], "#nu_{#mu} CC", "l")

    if len(cutNames)>1:
      lg_1.AddEntry(hs_Truth[1], "Cut 1. Fiducial volume", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>2:
      lg_1.AddEntry(hs_Truth[2], "Cut 2. Not a clear cosmic", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>3:
      lg_1.AddEntry(hs_Truth[3], "Cut 3. Flash matching score<12", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>4:
      lg_1.AddEntry(hs_Truth[4], "Cut 4. Side CRT veto", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>5:
      lg_1.AddEntry(hs_Truth[5], "Cut 5. Top CRT veto", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>6:
      lg_1.AddEntry(hs_Truth[6], "Cut 6. Muon track found", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    if len(cutNames)>7:
      lg_1.AddEntry(hs_Truth[7], "Cut 6. Contained muon track found", "l")
    else:
      lg_1.AddEntry(h_empty_legend, "", "")

    lg_1.Draw()

    ## ratio

    c1_down.cd()

    h_dummy_down = ROOT.TH1D('h_dummy_down', '', int( (xMax-xMin)/dx ), xMin, xMax)
    h_dummy_down.GetXaxis().SetTitle(xTitle)
    h_dummy_down.Draw("axis")
    h_dummy_down.GetYaxis().SetTitle("Efficiency")

    dict_h_effs = dict()
    list_h_eff_Total = []
    for i_h in range(1,len(hs_Truth)):
      h_eff = hs_Truth[i_h].Clone('hs_Truth_eff_%d_over_%d'%(i_h,i_h-1))
      h_eff_Total = hs_Truth[i_h].Clone('hs_Truth_eff_%d_Total'%(i_h))
      dict_h_effs['hs_Truth_eff_%d_over_%d'%(i_h,i_h-1)] = h_eff
      dict_h_effs['hs_Truth_eff_%d_Total'%(i_h)] = h_eff_Total

      if i_h!=0:
        print( '%f\t%f' % ( h_eff.Integral()/hs_Truth[0].Integral(), h_eff.Integral()/hs_Truth[i_h-1].Integral() ) )

      h_eff.Divide(hs_Truth[i_h-1])
      if i_h==len(hs_Truth)-1:
        h_eff.Draw("histsame")
      #h_eff.Draw("histsame")

      h_eff_Total.Divide(hs_Truth[0])
      list_h_eff_Total.append(h_eff_Total)

    ## Finalize

    h_dummy_up, h_dummy_down = canvas_margin.hist_axis(h_dummy_up, h_dummy_down)

    c1.cd()

    latex_POT.DrawLatex(0.70, 0.96, 'NuMI, '+str_TargetPOT);

    if varName=="MuonTruthNuMICosineTheta":
      lt.DrawLatex(0.30, 0.45, "#splitline{ICARUS Simulation}{Preliminary}")
    else:
      lt.DrawLatex(0.60, 0.45, "#splitline{ICARUS Simulation}{Preliminary}")

    c1.SaveAs(outDirBase+'/'+sample+'_'+varName+'_TrueEfficiency.pdf')

    ## Efficiency plot

    c1_up.cd()
    h_dummy_up.Draw('hist')
    h_dummy_up.GetYaxis().SetRangeUser(0.,1.5)
    h_dummy_up.GetYaxis().SetTitle("Cumulative efficiency")
    for h in list_h_eff_Total:
      h.Draw('histsame')
    lg_1.Draw()

    c1.SaveAs(outDirBase+'/'+sample+'_'+varName+'_TrueEfficiency_Total.pdf')
    c1.Close()

