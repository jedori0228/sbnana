#pragma once

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
#include "TGraphAsymmErrors.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

#include "TObject.h"
#include <iostream>

//==== Systematic
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

using namespace std;

namespace ICARUSNumuXsec{

  class HistoProducer{

  public:

    //==== Default Constructor with p4
    HistoProducer();

    void initialize();
    bool setCut(TString cutName);
    void bookSlice(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
    void bookTruth(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
    void bookRecoMuon(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
    void bookRecoProton(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);
    void bookRecoNeutrino(SpectrumLoader& loader, Cut cut, SpillCut spillCut=kNoSpillCut);

    void saveHistograms();

    //==== Systematic weights
    std::vector<const ISyst*> IGENIESysts;
    std::vector<Var> systGENIEWs;
    void setSystematicWeights();
    std::vector<const ISyst*> IFluxSysts;

    ~HistoProducer();

    double TargetPOT;
    TString str_TargetPOT;

    TString outputDir;
    TFile *outputfile;
    TString outputName;

    TString currentCutName;
    vector<TString> vec_cutNames;
    map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;
    //==== EnsembleSpectrum-based systematics
    map< TString, vector<EnsembleSpectrum *> > map_cutName_to_vec_EnsembleSpectrums;
    //==== Other systematics
    map< TString, vector< pair<TString, Spectrum *> > > map_cutName_to_vec_SystSpectrumPairs;

    //==== booleans
    bool doTruthMatch;
    bool doPerformanceStudy;

  private:

  };

}

