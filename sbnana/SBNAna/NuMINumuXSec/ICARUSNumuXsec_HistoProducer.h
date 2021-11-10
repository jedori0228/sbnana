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
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Systematics.h"

using namespace std;

namespace ICARUSNumuXsec{

  class HistoProducer{

  public:

    //==== Default Constructor with p4
    HistoProducer();

    void initialize();
    bool setCut(TString cutName);
    void bookSpectrums(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookSlice(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookTruth(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookRecoMuon(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookMuonPerformance(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookRecoProton(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void bookRecoNeutrino(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);

    void saveHistograms();

    //==== Systematic weights
    std::vector<const ISyst*> IGENIESysts;
    std::vector< vector<Var> > vec_UniverseWeightsForEachGENIESource;
    void setSystematicWeights();
    std::vector<const ISyst*> IFluxSysts;
    std::vector<const ISyst*> IDetectorSysts;

    ~HistoProducer();

    double TargetPOT;
    TString str_TargetPOT;

    TString outputDir;
    TFile *outputfile;
    TString outputName;

    TString currentCutName;
    vector<TString> vec_cutNames;
    std::map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;
    //==== EnsembleSpectrum-based systematics
    std::map< TString, vector< pair<TString, EnsembleSpectrum *> > > map_cutName_to_vec_SystEnsembleSpectrumPairs;
    void addEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);
    void addEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);

    //==== Other systematics
    std::map< TString, vector< pair<TString, Spectrum *> > > map_cutName_to_vec_SystSpectrumPairs;
    void addUpDownSystematic(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);
    void addUpDownSystematic(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);

    //==== booleans
    bool doTruthMatch;
    bool doPerformanceStudy;
    bool fillBeamInfo;
    bool fillNominal;

  private:

  };

}

