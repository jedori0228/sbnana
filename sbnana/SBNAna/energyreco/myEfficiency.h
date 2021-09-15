#ifndef myEfficiency_h
#define myEfficiency_h

// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include "sbnana/SBNAna/energyreco/myMuonSelection.h"
#include "sbnana/SBNAna/energyreco/myProtonSelection.h"
#include "sbnana/SBNAna/energyreco/myEstimator.h"
#include "sbnana/SBNAna/energyreco/myTruth.h"
#include "sbnana/SBNAna/energyreco/myEventSelection.h"

#include "sbnana/SBNAna/energyreco/myFilelist.h"
#include "sbnana/SBNAna/energyreco/myFilelistShort.h"
#include "sbnana/SBNAna/energyreco/myFilelistReCAF.h"

#include "TObject.h"
#include "TVector3.h"
#include <iostream>

//==== Systematic
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

using namespace std;

class myEfficiency{

public:

  //==== Default Constructor with p4
  myEfficiency();

  void initialize();
  void bookTruth(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
  void bookRecoMuon(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
  void bookRecoProton(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
  void bookRecoNeutrino(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);

  void saveHistograms();

  //==== Systematic weights
  std::vector<const ISyst*> ISysts;
  std::vector<Var> systWs;
  void setSystematicWeights();

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
