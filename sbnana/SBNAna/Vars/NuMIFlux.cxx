#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana {

  //// ----------------------------------------------
  NuMIPpfxFluxWeight::NuMIPpfxFluxWeight()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    fFluxFilePath = std::string(sbndata) +
                   "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";

    TFile f(fFluxFilePath.c_str());
    if (f.IsZombie()) {
      std::cout << "NuMIPpfxFluxWeight: Failed to open " << fFluxFilePath << std::endl;
      std::abort();
    }


    fFluxFilePath_G3Chase = std::string(sbndata) +
                   "beamData/NuMIdata/g3Chase_weights_rewritten.root";
    TFile f_G3Chase(fFluxFilePath_G3Chase.c_str());
    if (f_G3Chase.IsZombie()) {
      std::cout << "NuMIPpfxFluxWeight: Failed to open " << fFluxFilePath_G3Chase << std::endl;
      std::abort();
    }


    for (int hcIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          std::string hNamePPFX = "ppfx_flux_weights/hweights_";
          if (hcIdx == 0)
            hNamePPFX += "fhc_";
          else
            hNamePPFX += "rhc_";
          if (flavIdx == 0)
            hNamePPFX += "nue";
          else
            hNamePPFX += "numu";
          if (signIdx == 1) hNamePPFX += "bar";

          TH1* h_ppfx = (TH1*)f.Get(hNamePPFX.c_str());
          if (!h_ppfx) {
            std::cout << "NuMIPpfxFluxWeight: failed to find " << hNamePPFX << " in " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;

          // g3 chase
          std::string hNameG3Chase = "hweights_";
          if (hcIdx == 0)
            hNameG3Chase += "fhc_";
          else
            //hNameG3Chase += "rhc_"; //TODO
            hNameG3Chase += "fhc_"; //TODO
          if (flavIdx == 0)
            hNameG3Chase += "nue";
          else
            hNameG3Chase += "numu";
          if (signIdx == 1) hNameG3Chase += "bar";

          hNameG3Chase += "_total";

          TH1* h_G3Chase = (TH1*)f_G3Chase.Get(hNameG3Chase.c_str());
          if (!h_G3Chase) {
            std::cout << "NuMIPpfxFluxWeight: failed to find " << hNameG3Chase << " in " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_G3Chase = (TH1*)h_G3Chase->Clone(UniqueName().c_str());
          h_G3Chase->SetDirectory(0);

          fWeight_G3Chase[hcIdx][flavIdx][signIdx] = h_G3Chase;


        }
      }
    }
  }

  NuMIPpfxFluxWeight::~NuMIPpfxFluxWeight()
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k){
          delete fWeight[i][j][k];
          delete fWeight_G3Chase[i][j][k];
        }
  }

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeight([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxWeight(&slc->truth);
  });

  //// ----------------------------------------------

  const TruthVar kGetTruthNuMIFluxWeight([](const caf::SRTrueInteractionProxy* nu) -> double {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    if (!FluxWeightNuMI.fWeight[0][0][0]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;

    TH1* h = FluxWeightNuMI.fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(nu->E);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return 1.0;
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) return 1.0;
    return h->GetBinContent(bin);
  });

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeightWithG3Chase([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxWeightWithG3Chase(&slc->truth);
  });

  //// ----------------------------------------------

  const TruthVar kGetTruthNuMIFluxWeightWithG3Chase([](const caf::SRTrueInteractionProxy* nu) -> double {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    if (!FluxWeightNuMI.fWeight[0][0][0]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;

    TH1* h = FluxWeightNuMI.fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(nu->E);

    double w_ppfx = 1.;

    if (bin == 0 || bin == h->GetNbinsX() + 1) w_ppfx = 1.0;
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) w_ppfx = 1.0;
    w_ppfx = h->GetBinContent(bin);

    // G3Chase

    TH1* h_G3Chase = FluxWeightNuMI.fWeight_G3Chase[hcIdx][flavIdx][signIdx];
    assert(h_G3Chase);

    const int bin_G3Chase = h_G3Chase->FindBin(nu->E);

    double w_G3Chase = 1.;

    if (bin_G3Chase == 0 || bin_G3Chase == h_G3Chase->GetNbinsX() + 1) w_G3Chase = 1.0;
    if (std::isinf(h_G3Chase->GetBinContent(bin_G3Chase)) || std::isnan(h_G3Chase->GetBinContent(bin_G3Chase))) w_G3Chase = 1.0;
    w_G3Chase = h_G3Chase->GetBinContent(bin_G3Chase);


    return w_ppfx * w_G3Chase;





  });

}
