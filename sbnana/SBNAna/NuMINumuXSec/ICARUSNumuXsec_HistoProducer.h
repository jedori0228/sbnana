#pragma once

#include <iostream>

// ROOT
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

// NuMINumuXSec
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Weights.h"

// ThreeTrack
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_ThreeTrack_Cuts.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Cuts.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_MichelStudy.h"

using namespace ana;
using namespace std;
using namespace ICARUSNumuXsec;

namespace ICARUSNumuXsec{

  class HistoProducer{

  public:

    HistoProducer();

    void initialize();
    bool setCut(TString cutName);

    // Spectrum booking
    // - helpers
    void FillFlashMatching(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void FillLongestTrack(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void FillSpill(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void FillSlice(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 221121_CRTPMTMatching
    void CRTPMTMatchingStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - NuMu event selection
    void NuMINumuXSec(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230216_ThreeTrack
    void TruthStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void TrackPIDStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void HighChi2MuonProton(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void TwoTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void TwoTrackTruthMatching(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void ThreeTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230418_StubStudy
    void StubStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230517_TriggerEffStudy
    void TriggerEffStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230524_MichelStudy
    void MichelStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);

    void Test(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);

    void saveHistograms();

    // Systematic weights
    std::string SystProviderPrefix;
    void setSystematicWeights();
    std::vector<const ISyst*> IAllSysts; // TODO update later
    std::vector<const ISyst*> IGENIESysts;

    std::map<std::string, std::vector<Var>> map_DepDialName_to_UniverseWeights;
    std::map<std::string, std::vector<SpillVar>> map_DepDialName_to_UniverseSpillWeights;

    int NNuMIFluxPCA;
    std::vector<const ISyst*> IFluxSysts;

    std::vector<const ISyst*> IDetectorSysts;

    ~HistoProducer();

    double TargetPOT;
    TString str_TargetPOT;

    TString outputDir;
    TFile *outputfile;
    TString outputName;
    TString sampleName;

    TString currentCutName;
    vector<TString> vec_cutNames;
    std::map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;

    // EnsembleSpectrum-based systematics
    std::map< TString, vector< pair<TString, EnsembleSpectrum *> > > map_cutName_to_vec_SystEnsembleSpectrumPairs;
    void AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);
    void AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);

    // Systematics by ISyst
    std::map< TString, vector< pair<TString, Spectrum *> > > map_cutName_to_vec_SystSpectrumPairs;
    void AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);
    void AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);

    // slice
    template<class T>
    void FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& var, const Binning& binning, SpillCut spillCut, Cut cut, bool ForceFill=false);
    template<class T>
    void FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& var1, const Binning& binning1, const T& var2, const Binning& binning2, SpillCut spillCut, Cut cut, bool ForceFill=false);
    // spill
    template<class T>
    void FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& spillvar, const Binning& binning, SpillCut spillCut, bool ForceFill=false);
    template<class T>
    void FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& spillvar1, const Binning& binning1, const T& spillvar2, const Binning& binning2, SpillCut spillCut, bool ForceFill=false);

    void FillSystSpectrum(SpectrumLoader& loader, const std::string& label, const Var& var, const Binning& binning, SpillCut spillCut, Cut cut);
    void FillSystSpectrum(SpectrumLoader& loader, const std::string& label, const SpillVar& var, const Binning& binning, SpillCut spillCut);
    void FillCVandSystSpectrum(SpectrumLoader& loader, const std::string& label, const Var& var, const Binning& binning, SpillCut spillCut, Cut cut);
    void FillCVandSystSpectrum(SpectrumLoader& loader, const std::string& label, const SpillVar& var, const Binning& binning, SpillCut spillCut);

    // weights
    bool ApplyNuMIPPFXCVWeight;
    const Var GetGlobalWeight();

    // booleans
    bool IsData;
    bool FillMetaData;
    bool FillSystematics;
    bool FillGENIESyst;
    bool FillFlux;


  private:

  };

}

