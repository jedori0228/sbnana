#pragma once

#include "cetlib/search_path.h"

#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TH1.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TProfile.h"
#include "TFile.h"

#include <cassert>
#include <cstdlib>

#include <iostream>

#define M_MUON 0.1057
#define M_CHARGEDPION 0.13957039
#define M_NEUTRALPION 0.1349768
#define M_PIZERO 0.1349768
#define M_PROTON 0.938272
#define M_NEUTRON 0.939565
#define M_ELECTRON 0.00051
#define E_EffNuclB 0.040

using namespace std;

namespace ICARUSNumuXsec{

  //---------------------------------------------------
  class VolumeTool{

  public:

    VolumeTool()
    {
      fvCryo0 = {-358.49, -61.94,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.95};
      fvCryo1 = { +61.94, +358.49,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.95};
    }

    FidVol fvCryo0, fvCryo1;
    mutable double XMargin, YMargin, ZMarginUp, ZMarginDown;

    bool isContained(double x, double y, double z) const;
    int containedCryo(double x, double y, double z) const;

    int TPCIndex(double x, double y, double z) const;

    double GetTotalVolume() const;


  }; // END Class VolumeTool

  //---------------------------------------------------
  class ActiveVolumeTool : public VolumeTool{

  public:

    ActiveVolumeTool(){
      XMargin = 0.;
      YMargin = 0.;
      ZMarginUp = 0.;
      ZMarginDown = 0.;
    }

    static ActiveVolumeTool& Instance();

  }; // END ActiveVolumeTool class def

  //---------------------------------------------------
  class VertexContained : public VolumeTool{

  public:

    VertexContained(){
      XMargin = 25.;
      YMargin = 25.;
      ZMarginUp = 30.;
      ZMarginDown = 50.;
    }

    static VertexContained& Instance();

  }; // END VertexContained class def

  //---------------------------------------------------
  class TrackContained : public VolumeTool{

  public:

    TrackContained(){
      XMargin = 10.;
      YMargin = 10.;
      ZMarginUp = 10.;
      ZMarginDown = 10.;
    }

    static TrackContained& Instance();

  }; // END TrackContained class def

  //---------------------------------------------------
  class NuMICoordinateTool{

  public:

    NuMICoordinateTool();

    TVector3 GetICARUSCoord(double x, double y, double z) const;

    static NuMICoordinateTool& Instance();

  protected:

    TVector3 NuDirection_NuMI;
    TMatrixD rotMatNtoI = TMatrixD(3,3);
    TMatrixD tranVecNtoI = TMatrixD(3,1);

  }; // END Class NuMICoordinateTool

  //---------------------------------------------------
  class dEdXTemplateTool{

  public:

    dEdXTemplateTool();

    double GetdEdX(double rr, int ptlType) const;
    double GetdEdXErr(double rr, int ptlType) const;

    double CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo, int ptlType) const;
    double CalculateInelasticPionChi2(const caf::Proxy<caf::SRTrackCalo>& calo, int ptlType) const;

    static dEdXTemplateTool& Instance();

    std::string fTemplateFile;
    std::string fROOTfile;

    TProfile *dedx_range_pro;   ///< proton template
    TProfile *dedx_range_ka;    ///< kaon template
    TProfile *dedx_range_pi;    ///< pion template
    TProfile *dedx_range_mu;    ///< muon template

  };

  //---------------------------------------------------
  class SterileNuTool{

  public:

    SterileNuTool();

    double GetOscProb(double LoE, int i, int f) const;

    static SterileNuTool& Instance();

    double sin2th;
    double m2;

  };

  //---------------------------------------------------
  class ParticleTool{

  public:

    ParticleTool();

    static ParticleTool& Instance();
    double GetMass(int pdg) const;

  };

  //---------------------------------------------------
  class InteractionTool{

    public:

      struct NParticles{
        int NMuon = 0;
        int NProton = 0;
        int NNeutron = 0;
        int NPip = 0;
        int NPim = 0;
        int NPi0 = 0;
      };

      InteractionTool();

      static InteractionTool& Instance();

      mutable bool UseGHepRecord;

      mutable vector<int> MuonIndices;
      mutable vector<int> ProtonIndices;
      mutable vector<int> NeutronIndices;
      mutable vector<int> PipIndices;
      mutable vector<int> PimIndices;
      mutable vector<int> Pi0Indices;
      void ClearIndices() const;

      bool MuonContained() const;

      NParticles GetNParticles(const caf::SRSliceProxy* slc) const;
      NParticles GetNParticles(const caf::Proxy<std::vector<caf::SRTrueParticle>>& prim) const;

  };

  //---------------------------------------------------
  // Printing
  void PrintPrimaries(const caf::SRSliceProxy* slc);
  // For a given truth particle, find the reco object
  int GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx, double scorecut=0.5);
  int GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx, double scorecut=0.5);
  std::vector<double> GetMatchedRecoShowerIndices(const caf::SRSliceProxy* slc, int truth_idx, double scorecut=0.5);
  int GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx);

  //---------------------------------------------------
  std::vector<std::string> GetGENIEMultisigmaKnobNames();
  std::vector<std::string> GetGENIEMorphKnobNames();
  std::vector<std::string> GetGENIEDependentKnobNames();
  std::vector<std::string> GetGENIEMultisimKnobNames();

  //---------------------------------------------------
  bool IsPFPTrack(const caf::SRPFPProxy& pfp);
  bool IsPFPShower(const caf::SRPFPProxy& pfp);

}
