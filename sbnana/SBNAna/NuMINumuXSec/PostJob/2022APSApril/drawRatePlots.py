import os,ROOT
from array import array
import tdrstyle
import canvas_margin
from mylib import *
import InTimeCosmicSpillCalculator

tdrstyle.setTDRStyle()
ROOT.TH1.AddDirectory(False)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetTitleXOffset(0.85)

## Normalization in this plot
POTtoNormalize = 6e20
str_POTtoNormalize = "6e20 POT"

print("@@ POTtoNormalize = %1.2e"%(POTtoNormalize))

latex_POT = ROOT.TLatex()
latex_POT.SetNDC()
latex_POT.SetTextSize(0.035)

## inputs

inputBaseDir = "/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/HistoProducer/"

inputDirName_Data = "Run7568_NUMI/v09_45_00/"
inputDirName_BeamMC = 'NUMI_Nu_Cosmics/v09_45_00/'
inputDirName_IntimeCosmicMC = 'NUMI_in-time_Cosmics/v09_45_00/'

inputDirName = "FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

commonFname = "output__Nominal__Incl_PerInt__CutFlow.root"

f_Data = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_Data+'/'+commonFname)
f_BeamMC = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_BeamMC+'/'+commonFname)
f_IntimeCosmicMC = ROOT.TFile(inputBaseDir+inputDirName+'/'+inputDirName_IntimeCosmicMC+'/'+commonFname)

## Data info
##   This is necessary to undo the scaling done by ToTH1
##   https://github.com/SBNSoftware/sbnana/blob/b085dfe7238dd31ce4ff7b2e791db5f313ed5ac5/sbnana/CAFAna/Core/Spectrum.cxx#L442
h_Data_TargetPOT = f_Data.Get("hist_TargetPOT")
Data_TargetPOT = h_Data_TargetPOT.GetBinContent(1)
Data_RealPOTFromTable = 1644717.858e12
nominalIntensity = 48.09116545e12
#nominalIntensity = 6e13

## For the in-time cosmic MC normalization
m_InTimeCosmicSpillCalculator = InTimeCosmicSpillCalculator.InTimeCosmicSpillCalculator(POTtoNormalize, nominalIntensity)
InTimeCosmicSpill = m_InTimeCosmicSpillCalculator.GetInTimeCosmicSpill(f_BeamMC, f_IntimeCosmicMC)
## https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=20331&filename=20201215-Guglielmi-For_today_discussion.pdf&version=1

## outfile
outDirBase = './plots/'+inputDirName+'/RatePlots/'+inputDirName_Data+'/'
os.system('mkdir -p '+outDirBase)

cuts = [

['NoCut',
 'No cut',
  [
'NeutrinoTruthE',
'FMScore',
'FMTime',
'NuScore',
'SliceTrackNhitsPlane0',
'SliceTrackNhitsPlane1',
'SliceTrackNhitsPlane2',
'SliceTrackChargePlane0',
'SliceTrackChargePlane1',
'SliceTrackChargePlane2',
'AllTrackStartPositionX',
'AllTrackStartPositionY',
'AllTrackStartPositionZ',
  ],
  [

  ],
],

['NotClearCosmic',
 'Not a clear cosmic',
  [
'FMScore'
  ],
  [

  ],
],

['HasMuon',
 'Muon track found',
  [
'NeutrinoTruthE',
'VertexRecoX',
'VertexRecoY',
'VertexRecoZ',
'FMScore',
'FMTime',
'NuScore',
'MuonLength',
'MuonRecoP',
'MuonRecoCosineTheta',
'MuonRecoNuMICosineTheta',
  ],
  [

  ],
],

['HasMuon_Contained',
 'Contained muon track found',
  [
'NeutrinoTruthE',
'VertexRecoX',
'VertexRecoY',
'VertexRecoZ',
'FMScore',
'FMTime',
'NuScore',
'MuonLength',
'MuonRecoP',
'MuonRecoCosineTheta',
'MuonRecoNuMICosineTheta',
  ],
  [

  ],
],

['HasMuon_Exiting',
 'Exiting muon track found',
  [
'NeutrinoTruthE',
'VertexRecoX',
'VertexRecoY',
'VertexRecoZ',
'FMScore',
'FMTime',
'NuScore',
'MuonLength',
'MuonRecoP',
'MuonRecoCosineTheta',
'MuonRecoNuMICosineTheta',
  ],
  [

  ],
],

['NoCRT_HasMuon',
 'Muon track found (No CRT)',
  [
'NeutrinoTruthE',
'VertexRecoX',
'VertexRecoY',
'VertexRecoZ',
'FMScore',
'FMTime',
'NuScore',
'MuonRecoP',
'MuonRecoCosineTheta',
'MuonRecoNuMICosineTheta',
'MuonLength',
'FMScore',
'FMTime',
'NuScore',
  ],
  [

  ],
],

['NoCRT_HasMuon_Contained',
 'Contained muon track found (No CRT)',
  [
'NeutrinoTruthE',
'VertexRecoX',
'VertexRecoY',
'VertexRecoZ',
'FMScore',
'FMTime',
'NuScore',
'MuonRecoP',
'MuonRecoCosineTheta',
'MuonRecoNuMICosineTheta',
'MuonLength',
'NeutrinoCombinedEnergy',
'FMScore',
'FMTime',
'NuScore',
'ProtonRecoP',
#'ProtonLength',
  ],
  [

  ],
],

]

samples = [
[f_IntimeCosmicMC, 'AllSamples', 'In-time cosmic', ROOT.kOrange, "bkgd"],

#[f_BeamMC, 'Cosmic', 'Out-time cosmic', ROOT.kYellow+2, "bkgd"],
[f_BeamMC, 'Cosmic', 'Cosmic background', ROOT.kOrange, "bkgd"],


#[f_BeamMC, 'NuECC', '#nu_{e}-CC', ROOT.kGreen, "bkgd"],

[f_BeamMC, 'NuMuNC', '#nu_{#mu}-NC', ROOT.kBlue, "bkgd"],

#[f_BeamMC, 'NuMuCC', '#nu_{#mu}-CC', ROOT.kRed, "signal"],
[f_BeamMC, 'NuMuCCSignalDef', '#nu_{#mu}-CC', ROOT.kRed, "signal"],

]

vec_P_bins = [0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
vec_CosineTheta_bins = [-1, -0.7, -0.5, -0.2, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

## y=1 graph
g1_x = [-9000, 9000]
g1_y = [1, 1]
g1 = ROOT.TGraph(2, array("d", g1_x ), array("d", g1_y ))

for cut in cuts:

  cutName = cut[0]
  cutAlias = cut[1]
  varNames_1D = cut[2]
  varNames_2D = cut[3]

  outDir = outDirBase+'/'+cutName+'/'
  os.system('mkdir -p '+outDir)
  f_out = ROOT.TFile(outDir+"/hists_RecoVariable.root","RECREATE")

  ## 1D

  for varName in varNames_1D:

    ## Efficiency based on true disbtribution

    xTitle = ''
    xMin = 0.
    xMax = 5.
    dx = 0.1
    nRebin = 1

    if varName=='MuonRecoP':
      xTitle = 'Muon momentum (GeV)'
      xMin = 0.
      xMax = 3.
      dx = 0.1
      nRebin = 2
    elif varName=='MuonRecoCosineTheta':
      xTitle = 'Muon cos(#theta)_{BNB}'
      xMin = -1.
      xMax = 1.
      dx = 0.1
    elif varName=='MuonRecoNuMICosineTheta':
      xTitle = 'Muon cos(#theta_{NuMI})'
      xMin = -1.
      xMax = 1.
      dx = 0.1
    elif varName=='MuonLength':
      xTitle = 'Muon length (cm)'
      xMin = 0.
      xMax = 600.
      dx = 1
      nRebin = 50
    elif varName=='ProtonRecoP':
      xTitle = 'Proton momentum (GeV)'
      xMin = 0.
      xMax = 5.
      dx = 0.1
      nRebin = 1
    elif varName=='ProtonLength':
      xTitle = 'Proton length (cm)'
      xMin = 0.
      xMax = 200.
      dx = 1.
      nRebin = 10
    elif varName=='FMScore':
      xTitle = 'Flash matching score'
      xMin = -2.
      xMax = 100.
      dx = 1
      nRebin = 1
      if cutName=='NoCut':
        xMin = -1.
        xMax = 100.
      elif cutName=="NotClearCosmic":
        xMin = 0.
        xMax = 40.
    elif varName=='FMTime':
      xTitle = 'Flash matching time'
      xMin = -50.
      xMax = 50.
      dx = 1
      nRebin = 1
    elif varName=='NuScore':
      xTitle = '#nu-score'
      xMin = 0.
      xMax = 1.
      dx = 0.01
      nRebin = 1
    elif varName=='NeutrinoCombinedEnergy':
      xTitle = 'Neutrino energy (GeV)'
      xMin = 0.
      xMax = 10.
      dx = 0.1
      nRebin = 5
    elif varName=='SliceTrackNhitsPlane0':
      xTitle = 'Number of thits of tracks (Ind. 1)'
      xMin = 0.
      xMax = 2000.
      dx = 1
      nRebin = 1
    elif varName=='SliceTrackNhitsPlane1':
      xTitle = 'Number of hits of tracks (Ind. 2)'
      xMin = 0.
      xMax = 2000.
      dx = 1
      nRebin = 1
    elif varName=='SliceTrackNhitsPlane2':
      xTitle = 'Number of hits of tracks (Coll.)'
      xMin = 0.
      xMax = 2000.
      dx = 1
      nRebin = 1
    elif varName=='SliceTrackChargePlane0':
      xTitle = 'Charge of tracks (Ind. 1), #times 10^{3}'
      xMin = 0.
      xMax = 10000.
      dx = 1
      nRebin = 10
    elif varName=='SliceTrackChargePlane1':
      xTitle = 'Charge of tracks (Ind. 2), #times 10^{3}'
      xMin = 0.
      xMax = 10000.
      dx = 1
      nRebin = 10
    elif varName=='SliceTrackChargePlane2':
      xTitle = 'Charge of tracks (Coll.), #times 10^{3}'
      xMin = 0.
      xMax = 10000.
      dx = 1
      nRebin = 10
    elif varName=='AllTrackStartPositionX':
      xTitle = 'Track start position, x (cm)'
      xMin = -500
      xMax = 500.
      dx = 1
      nRebin = 1
    elif varName=='AllTrackStartPositionY':
      xTitle = 'Track start position, y (cm)'
      xMin = -200
      xMax = 200.
      dx = 1
      nRebin = 1
    elif varName=='AllTrackStartPositionZ':
      xTitle = 'Track start position, z (cm)'
      xMin = -1000
      xMax = 1000.
      dx = 1
      nRebin = 1
    elif varName=='FMChargeCenterX':
      xTitle = 'FlashMatching charge center, x (cm)'
      xMin = -500
      xMax = 500.
      dx = 1
      nRebin = 1
    elif varName=='FMChargeCenterY':
      xTitle = 'FlashMatching charge center, y (cm)'
      xMin = -200
      xMax = 200.
      dx = 1
      nRebin = 1
    elif varName=='FMChargeCenterZ':
      xTitle = 'FlashMatching charge center, z (cm)'
      xMin = -1000
      xMax = 1000.
      dx = 1
      nRebin = 1
    elif varName=='FMLightCenterX':
      xTitle = 'FlashMatching light center, x (cm)'
      xMin = -500
      xMax = 500.
      dx = 1
      nRebin = 1
    elif varName=='FMLightCenterY':
      xTitle = 'FlashMatching light center, y (cm)'
      xMin = -200
      xMax = 200.
      dx = 1
      nRebin = 1
    elif varName=='FMLightCenterZ':
      xTitle = 'FlashMatching light center, z (cm)'
      xMin = -1000
      xMax = 1000.
      dx = 1
      nRebin = 1
    elif varName=='FMChargeQ':
      xTitle = 'FlashMatching charge, #times 10^{3}'
      xMin = -2
      xMax = 500
      dx = 1
      nRebin = 1
    elif varName=='MuonChi2Muon':
      xTitle = '#chi^{2}_{#mu} of muon track'
      xMin = 0
      xMax = 150
      dx = 1
      nRebin = 10
    elif varName=='MuonChi2Proton':
      xTitle = '#chi^{2}_{proton} of muon track'
      xMin = 0
      xMax = 150
      dx = 1
      nRebin = 10
    elif varName=='ProtonChi2Proton':
      xTitle = '#chi^{2}_{proton} of proton track'
      xMin = 0
      xMax = 150
      dx = 1
      nRebin = 10
    elif varName=='VertexRecoX':
      xTitle = 'Track start position, x (cm)'
      xMin = -500
      xMax = 500.
      dx = 1
      nRebin = 20
    elif varName=='VertexRecoY':
      xTitle = 'Track start position, y (cm)'
      xMin = -200
      xMax = 200.
      dx = 1
      nRebin = 10
    elif varName=='VertexRecoZ':
      xTitle = 'Track start position, z (cm)'
      xMin = -1000
      xMax = 1000.
      dx = 1
      nRebin = 50

    ## Rate

    c_Rate = ROOT.TCanvas('c_Rate', '', 800, 800)
    c_Rate.cd()

    hist_dummy_Rate = ROOT.TH1D('hist_dummy_Rate', '', int( (xMax-xMin)/dx ), xMin, xMax)
    hist_dummy_Rate.Draw("axis")
    hist_dummy_Rate.GetXaxis().SetTitle(xTitle)

    lg_Rate = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    lg_Rate.SetBorderSize(0)
    lg_Rate.SetFillStyle(0)

    ## data
    h_Data = f_Data.Get('AllSamples_'+cutName+"/"+varName+'_AllSamples_'+cutName)
    h_Data.Rebin(nRebin)
    h_Data.SetMarkerStyle(20)
    h_Data.SetMarkerColor(ROOT.kBlack)
    ## InTime cosmic is already scaled by POTtoNormalize/1e18
    ## https://github.com/SBNSoftware/sbnana/blob/e030b862847d8e84dc35006906276ca58be396ca/sbnana/CAFAna/Core/Spectrum.cxx#L442
    DataScaleFactor = (POTtoNormalize/Data_RealPOTFromTable) / (Data_TargetPOT/1e18)
    h_Data.Scale( DataScaleFactor )
    gr_Data = convertToGraph(h_Data)
    gr_Data.SetLineWidth(2)
    gr_Data.SetMarkerStyle(20)
    print("@@ DataScaleFactor",DataScaleFactor)
    print("@@ Data integral",h_Data.Integral())

    #lg_Rate.AddEntry(gr_Data, "Data", "ep")

    hs = []
    y_max = h_Data.GetMaximum()
    y_max = -1
    h_CosmicSum = 0
    for sample in samples:

      f_histo = sample[0]
      sampleName = sample[1]
      sampleAlias = sample[2]
      sampleColor = sample[3]
      sampleType = sample[4]

      dirName = sampleName+'_'+cutName
      tDir = f_histo.Get( dirName )

      h = tDir.Get( varName+'_'+dirName )
      if not h:
        continue
      h.SetLineColor(sampleColor)
      h.Rebin(nRebin)
      h.SetName(sampleAlias+'\t'+sampleType)

      if sampleAlias=='In-time cosmic':
        h_Livetime = tDir.Get('Livetime_'+dirName )
        this_Livetime = h_Livetime.GetBinContent(1)
        print("In-time cosmic MC contains %1.2e spills"%(this_Livetime))
        ## InTime cosmic is already scaled by IntimeCosmicMC_TargetPOT/1e18
        ## https://github.com/SBNSoftware/sbnana/blob/e030b862847d8e84dc35006906276ca58be396ca/sbnana/CAFAna/Core/Spectrum.cxx#L442
        SpillScale = InTimeCosmicSpill/this_Livetime
        print("-> to match InTimeCosmicSpill (%1.2e), In-time MC should be scaled by %1.2e"%(InTimeCosmicSpill, SpillScale))
        alreadyScaled = m_InTimeCosmicSpillCalculator.IntimeCosmicMC_TargetPOT/1e18
        print("-> In-time MC has been already scaled by %1.2e/1e18 = %1.2e"%(m_InTimeCosmicSpillCalculator.IntimeCosmicMC_TargetPOT, alreadyScaled))
        h.Scale(SpillScale/alreadyScaled)
        print("In-time cosmic scaled by %1.2e"%(SpillScale/alreadyScaled))
      else:
        h.Scale(POTtoNormalize/m_InTimeCosmicSpillCalculator.BeamMC_TargetPOT)
        print("Beam MC scaled by %1.2e"%(POTtoNormalize/m_InTimeCosmicSpillCalculator.BeamMC_TargetPOT))

      if "cosmic" in sampleAlias.lower():
        if h_CosmicSum==0:
          h_CosmicSum = h.Clone()
        else:
          h_CosmicSum.Add(h)
          h_CosmicSum.Draw("histsame")
          y_max = max(y_max, h_CosmicSum.GetMaximum())
          lg_Rate.AddEntry(h_CosmicSum, "Cosmic background", 'l')
          hs.append(h_CosmicSum)
      else:
        h.Draw('histsame')
        y_max = max(y_max, h.GetMaximum())
        if 'NC' not in sampleAlias:
          lg_Rate.AddEntry(h, sampleAlias, 'l')
        hs.append(h)

      print(cutName,sampleName,h.Integral())

    #if h_Data.Integral()>0:
    #  h_Data.Draw("psame")

    lg_Rate.Draw()
    hist_dummy_Rate.GetYaxis().SetRangeUser(0., 1.2*y_max)

    if 'bnb' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.16, 0.96, 'BNB, '+str_POTtoNormalize);
    elif 'numi' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.16, 0.96, 'NuMI, '+str_POTtoNormalize);

    c_Rate.SaveAs(outDir+varName+'.pdf')

    ## Normalized

    y_max = -999.
    hs_Normed = []
    hist_dummy_Rate.Draw('axis')
    for h in hs:

      h_Normed = h.Clone(h.GetName()+"_Normed")
      if h_Normed.Integral()>0:
        h_Normed.Scale(1./h_Normed.Integral())

      if 'NC' not in h.GetName():
        h_Normed.Draw('histsame')
        y_max = max(y_max, h_Normed.GetMaximum())

      hs_Normed.append(h_Normed)
    hist_dummy_Rate.GetYaxis().SetRangeUser(0., 1.2*y_max)
    if varName=='NuScore':
      if cutName=='NoCut':
         hist_dummy_Rate.GetYaxis().SetRangeUser(0., 0.45*1.2)
    hist_dummy_Rate.GetYaxis().SetTitle("Normalized")

    h_Data_Normed = h_Data.Clone('h_Data_Normed')
    if h_Data.Integral()>0:
      h_Data_Normed.Scale(1./h_Data_Normed.Integral())
      gr_Data_Normed = convertToGraph(h_Data_Normed)
      gr_Data_Normed.SetLineWidth(2)
      gr_Data_Normed.SetMarkerStyle(20)
      #h_Data_Normed.Draw("psame")

    lg_Rate.Draw()
    c_Rate.SaveAs(outDir+varName+'_Normed.pdf')
    c_Rate.Close()

    ## Stack

    c_Stack = ROOT.TCanvas('c_Stack', '', 800, 800)

    c_Stack_up = ROOT.TPad("c_Stack_up", "", 0, 0.25, 1, 1)
    c_Stack_down = ROOT.TPad("c_Stack_down", "", 0, 0, 1, 0.25)

    c_Stack, c_Stack_up, c_Stack_down = canvas_margin.canvas_margin(c_Stack, c_Stack_up, c_Stack_down)
    c_Stack.cd()
    c_Stack.Draw()
    c_Stack_up.Draw()
    c_Stack_down.Draw()

    hist_dummy_Stack_up = ROOT.TH1D('hist_dummy_Stack_up', '', int( (xMax-xMin)/dx ), xMin, xMax)
    hist_dummy_Stack_down = ROOT.TH1D('hist_dummy_Stack_down', '', int( (xMax-xMin)/dx ), xMin, xMax)

    hist_dummy_Stack_up = ROOT.TH1D('hist_dummy_Stack_up', '', int( (xMax-xMin)/dx ), xMin, xMax)
    hist_dummy_Stack_down = ROOT.TH1D('hist_dummy_Stack_down', '', int( (xMax-xMin)/dx ), xMin, xMax)
    hist_dummy_Stack_up, hist_dummy_Stack_down = canvas_margin.hist_axis(hist_dummy_Stack_up, hist_dummy_Stack_down)

    c_Stack_up.cd()
    hist_dummy_Stack_up.Draw('axis')

    lg_Stack = ROOT.TLegend(0.55, 0.6, 0.9, 0.9)
    if 'CosineTheta' in varName:
      lg_Stack = ROOT.TLegend(0.25, 0.6, 0.6, 0.9)
    lg_Stack.SetBorderSize(0)
    lg_Stack.SetFillStyle(0)

    h_AllStack = ROOT.THStack('h_AllStack', '')
    h_AllSum = 0

    for h in hs:
      h.SetLineWidth(1)
      h.SetFillColor( h.GetLineColor() )

      h_AllStack.Add(h)
      if h_AllSum==0:
        h_AllSum = h.Clone('h_AllSum')
      else:
        h_AllSum.Add(h)

    for i_h in range(0,len(hs)):
      h = hs[len(hs)-1-i_h]
      h.SetLineWidth(1)
      h.SetFillColor( h.GetLineColor() )
      lg_Stack.AddEntry(h, h.GetName(), 'f')

    h_AllSum.SetFillColor(0)
    h_AllSum.SetLineColor(ROOT.kBlack)
    hist_dummy_Stack_up.GetYaxis().SetRangeUser(0., 1.1*max(h_AllSum.GetMaximum(),h_Data.GetMaximum()))

    '''
    if 'HasProton' in cutName:
      if varName=='MuonRecoP':
        hist_dummy_Stack_up.GetYaxis().SetRangeUser(0., 7000.)
      elif 'CosineTheta' in varName:
        hist_dummy_Stack_up.GetYaxis().SetRangeUser(0., 20000.)
    else:
      if varName=='MuonRecoP':
        hist_dummy_Stack_up.GetYaxis().SetRangeUser(0., 84000.)
      elif 'CosineTheta' in varName:
        hist_dummy_Stack_up.GetYaxis().SetRangeUser(0., 84000.)
    '''

    h_AllStack.Draw("histsame")

    gr_Data.Draw("psame")
    hist_dummy_Stack_up.Draw('axissame')

    c_Stack_down.cd()

    hist_dummy_Stack_down.Draw("axis")
    hist_dummy_Stack_down.GetYaxis().SetRangeUser(0., 2.0)
    hist_dummy_Stack_down.GetYaxis().SetTitle("Data/MC")
    hist_dummy_Stack_down.GetXaxis().SetTitle(xTitle)
    hist_dummy_Stack_down.SetNdivisions(504,"Y")
    h_Ratio = h_Data.Clone("h_Ratio")
    h_Ratio.Divide(h_AllSum)
    gr_Ratio = convertToGraph(h_Ratio)
    gr_Ratio.SetLineWidth(2)
    gr_Ratio.SetMarkerStyle(20)
    gr_Ratio.Draw("psame")
    g1.Draw("same")

    c_Stack.cd()
    lg_Stack.Draw()
    if 'bnb' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.16, 0.96, 'BNB, '+str_POTtoNormalize)
    elif 'numi' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.16, 0.96, 'NuMI, '+str_POTtoNormalize)
    if 'CosineTheta' in varName:
      latex_POT.DrawLatex(0.21, 0.91, cutAlias)
    else:
      latex_POT.DrawLatex(0.62, 0.91, cutAlias)

    c_Stack.SaveAs(outDir+varName+'_Stack.pdf')
    c_Stack.Close()

    ## Stack but normalized

    c_StackNormed = ROOT.TCanvas('c_StackNormed', '', 800, 800)

    c_StackNormed.cd()

    hist_dummy_StackNormed = ROOT.TH1D('hist_dummy_StackNormed', '', int( (xMax-xMin)/dx ), xMin, xMax)
    hist_dummy_StackNormed.GetXaxis().SetTitle(xTitle)
    hist_dummy_StackNormed.GetXaxis().SetTitleOffset(1.10)
    hist_dummy_StackNormed.Draw('axis')

    lg_StackNormed = ROOT.TLegend(0.50, 0.65, 0.9, 0.9)
    if 'CosineTheta' in varName:
      lg_StackNormed = ROOT.TLegend(0.20, 0.65, 0.6, 0.9)
    lg_StackNormed.SetBorderSize(0)
    lg_StackNormed.SetFillStyle(0)

    h_AllStackNormed = ROOT.THStack('h_AllStackNormed', '')
    h_AllSumNormed = 0
    h_BkgdSumNormed = 0
    MCNormScaleFactor = h_Data.Integral()/h_AllSum.Integral()
    hsNormed = []
    for h in hs:
      h.SetLineWidth(1)
      h.SetFillColor( h.GetLineColor() )
      h_Scale = h.Clone(h.GetName())
      #h_Scale.Scale(MCNormScaleFactor)

      h_AllStackNormed.Add(h_Scale)

      if h_AllSumNormed==0:
        h_AllSumNormed = h_Scale.Clone('h_AllSumNormed')
      else:
        h_AllSumNormed.Add(h_Scale)
      hsNormed.append(h_Scale)

      if h.GetName().split('\t')[1]=='bkgd':
        if h_BkgdSumNormed==0:
          h_BkgdSumNormed = h_Scale.Clone('h_BkgdSumNormed')
        else:
          h_BkgdSumNormed.Add(h_Scale)

    for i_h in range(0,len(hsNormed)):
      h = hsNormed[len(hs)-1-i_h]
      h.SetLineWidth(1)
      h.SetFillColor( h.GetLineColor() )
      if "In-time cosmic" not in h.GetName():
        lg_StackNormed.AddEntry(h, h.GetName(), 'f')

    h_AllStackNormed.Draw("histsame")
    hist_dummy_StackNormed.Draw('axissame')

    hist_dummy_StackNormed.GetYaxis().SetRangeUser(0., 1.1*h_AllSumNormed.GetMaximum())

    c_StackNormed.cd()
    lg_StackNormed.Draw()
    if 'bnb' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.70, 0.96, 'BNB, '+str_POTtoNormalize)
    elif 'numi' in inputDirName_BeamMC.lower():
      latex_POT.DrawLatex(0.70, 0.96, 'NuMI, '+str_POTtoNormalize)

    if 'CosineTheta' in varName:
      latex_POT.DrawLatex(0.22, 0.91, cutAlias)
    else:
      latex_POT.DrawLatex(0.52, 0.91, cutAlias)

    ## Write histograms
    f_out.cd()
    h_Data.SetName(cutName+'_'+varName+'_Data')
    h_Data.Write()
    h_AllSumNormed.SetName(cutName+'_'+varName+'_AllSum')
    h_AllSumNormed.Write()
    h_BkgdSumNormed.SetName(cutName+'_'+varName+'_BkgdSum')
    h_BkgdSumNormed.Write()

    c_StackNormed.SaveAs(outDir+varName+'_StackNormed.pdf')
    c_StackNormed.Close()

    ## Same as above, but only cosmic

    c_CosmicOnlyStackNormed = ROOT.TCanvas('c_CosmicOnlyStackNormed', '', 800, 800)

    c_CosmicOnlyStackNormed.cd()

    hist_dummy_StackNormed.Draw('axis')

    lg_CosmicOnlyStackNormed = ROOT.TLegend(0.50, 0.7, 0.9, 0.9)
    if 'CosineTheta' in varName:
      lg_CosmicOnlyStackNormed = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
    lg_CosmicOnlyStackNormed.SetBorderSize(0)
    lg_CosmicOnlyStackNormed.SetFillStyle(0)

    lg_CosmicOnlyStackNormed.AddEntry(h_Data, "Neutrino data (Run7568)", "pe")

    h_CosmicOnlyStackNormed = ROOT.THStack('h_CosmicOnlyStackNormed', '')
    MCNormScaleFactor = h_Data.Integral()/h_AllSum.Integral()
    print('@@ MCNormScaleFactor = ',MCNormScaleFactor)
    hsCosmicOnlyNormed = []
    for h in hs:
      h.SetLineWidth(1)
      h.SetFillColor( ROOT.kOrange )
      h.SetLineColor( ROOT.kOrange )
      h_Scale = h.Clone(h.GetName())
      h_Scale.Scale(MCNormScaleFactor)

      if "cosmic" in h.GetName():
        h_CosmicOnlyStackNormed.Add(h_Scale)
      hsCosmicOnlyNormed.append(h_Scale)

      if len(hsCosmicOnlyNormed)==1:
        lg_CosmicOnlyStackNormed.AddEntry(h_Scale, "Cosmic bacgkround", "f")


    h_CosmicOnlyStackNormed.Draw("histsame")
    hist_dummy_StackNormed.Draw('axissame')
    h_Data.Draw("pesame")

    hist_dummy_StackNormed.GetYaxis().SetRangeUser(0., 1.1*max(h_AllSumNormed.GetMaximum(),h_Data.GetMaximum()))

    c_CosmicOnlyStackNormed.cd()
    lg_CosmicOnlyStackNormed.Draw()
    latex_POT.DrawLatex(0.70, 0.96, 'NuMI, XXX POT')

    if 'CosineTheta' in varName:
      latex_POT.DrawLatex(0.21, 0.91, cutAlias)
    else:
      latex_POT.DrawLatex(0.50, 0.91, cutAlias)

    c_CosmicOnlyStackNormed.SaveAs(outDir+varName+'_CosmicOnlyStackNormed.pdf')
    c_CosmicOnlyStackNormed.Close()




