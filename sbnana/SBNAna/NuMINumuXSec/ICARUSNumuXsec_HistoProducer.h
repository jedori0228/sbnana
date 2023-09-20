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
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"
#include "sbnana/CAFAna/Systs/CalorimetrySysts.h"


// SBNANA
#include "sbnana/SBNAna/Vars/NuMIFlux.h"
#include "sbnana/SBNAna/Cuts/NuMIRelaxedVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecTreeHelper.h"
#include "sbnana/SBNAna/Cuts/NuMIDetSystStudy.h"

// NuMINumuXSec
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Weights.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"

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
    // - NuMu event selection
    void NuMIXSec(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void TruthStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void TwoTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - Track PID study (NuMI approval at 2023 July)
    void TrackPIDStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230418_StubStudy
    void StubStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230517_TriggerEffStudy
    void TriggerEffStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230524_MichelStudy
    void MichelStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230814_MakeTree
    bool TrueTreeFilled;
    void MakeTree(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void NuMIXSecBkgdStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230908_PIDStudy
    void MakePIDStudyTree(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 230918_DetSyst
    void MakeDetSystStudyTree(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);

    void Test(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);

    void saveHistograms();

    // Systematic weights
    std::string SystProviderPrefix;
    void setSystematicWeights();
    std::vector<const ISyst*> IAllSysts; // TODO update later

    std::vector<std::string> genieMultisigmaKnobNames;
    std::vector<const ISyst*> IGENIESysts;
    std::vector<const ISyst*> IGENIEMorphSysts;
    std::vector<std::string> genieDependentKnobNames;
    std::map<std::string, std::vector<Var>> map_DepDialName_to_UniverseWeights;
    std::map<std::string, std::vector<TruthVar>> map_DepDialName_to_TruthUniverseWeights;

    std::vector<std::string> geant4DependentKnobNames;


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
    std::map< TString, std::vector<ana::Tree *> > map_cutName_to_vec_Trees;
    bool MakeGUNDAMTree; // very spcial case for GUNDAM usage
    enum NSigmasSaveMode
    {
      kVector=0,
      kSpline=1,
      kGraph=2,
      kTClonesArrays=3,
    };
    NSigmasSaveMode nSigmasSaveMode;
    std::map< TString, std::vector<ana::NSigmasTree *> > map_cutName_to_vec_NSigmasTrees;
    std::map< TString, std::vector<ana::NUniversesTree *> > map_cutName_to_vec_NUniversesTrees;

    // EnsembleSpectrum-based systematics
    std::map< TString, vector< pair<TString, EnsembleSpectrum *> > > map_cutName_to_vec_SystEnsembleSpectrumPairs;

    // Systematics by ISyst
    std::map< TString, vector< pair<TString, Spectrum *> > > map_cutName_to_vec_SystSpectrumPairs;

    // Metadata
    void FillMetadata(SpectrumLoader& loader);

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

