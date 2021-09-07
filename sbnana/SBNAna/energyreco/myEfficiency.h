#ifndef myEfficiency_h
#define myEfficiency_h

// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include "canvas_margin.h"

#include "myMuonSelection.h"
#include "myProtonSelection.h"
#include "myEstimator.h"
#include "myTruth.h"
#include "myEventSelection.h"

#include "myFilelist.h"
#include "myFilelistShort.h"

#include "TObject.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

class myEfficiency{

public:

  //==== Default Constructor with p4
  myEfficiency();

  void initialize();
  void bookTruth(SpectrumLoader& loader, Cut cut);
  void bookRecoMuon(SpectrumLoader& loader, Cut cut);
  void bookRecoProton(SpectrumLoader& loader, Cut cut);
  void bookRecoNeutrino(SpectrumLoader& loader, Cut cut);

  void saveHistograms();

  ~myEfficiency();

  double TargetPOT;
  TString str_TargetPOT;

  TString outputDir;
  TFile *outputfile;
  TString outputName;

  vector<Spectrum *> vec_Spectrums;


private:

};

#endif
