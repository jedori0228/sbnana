#include "CAFAna/Core/OscCurve.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/func/IOscCalculator.h"

#include <cassert>
#include <iostream>
#include <map>

#include "TH1.h"

#include "CAFAna/Core/OscCalcSterileApprox.h"

namespace ana
{
  //----------------------------------------------------------------------
  OscCurve::OscCurve(osc::IOscCalculator* calc, int from, int to)
    : fFrom(from), fTo(to)
  {
    DontAddDirectory guard;

    fHist = HistCache::New(";True Energy (GeV);Probability", kTrueEnergyBins);

    // We have extra knowledge that calculators of this type have a special
    // mode allowing an intrinsic energy smearing.
    OscCalcSterileApprox* approx = DowncastToSterileApprox(calc, true);

    if(approx){
      for(int i = 0; i < fHist->GetNbinsX()+2; ++i){
        const double E = fHist->GetBinCenter(i);
        const double Elo = fHist->GetXaxis()->GetBinLowEdge(i);
        const double Ehi = fHist->GetXaxis()->GetBinUpEdge(i);
        // Use 2% resolution (intended to be << the resolution of any actual
        // event) or the bin width, whichever is larger
        fHist->SetBinContent(i, approx->P(from, to,
                                          std::min(Elo, 0.98*E),
                                          std::max(Ehi, 1.02*E)));
        fHist->SetBinError(i, 0);
      }
    }
    else{
      for(int i = 0; i < fHist->GetNbinsX()+2; ++i){
        const double E = fHist->GetBinCenter(i);
        fHist->SetBinContent(i, E > 0 ? calc->P(from, to, E) : 0);
        fHist->SetBinError(i, 0);
      }
    }
  }

  //----------------------------------------------------------------------
  OscCurve::OscCurve(TH1* h)
  {
    DontAddDirectory guard;

    const TString className = h->ClassName();

    if(className == "TH1D"){
      // Shortcut if types match
      fHist = HistCache::Copy((TH1D*)h);
    }
    else{
      fHist = HistCache::New("", h->GetXaxis());
      fHist->Add(h);
    }

    fHist->SetTitle(";True Energy (GeV);Probability");
  }

  //----------------------------------------------------------------------
  OscCurve::~OscCurve()
  {
    if(fHist && fHist->GetDirectory()){
      static bool once = true;
      if(once){
        once = false;
        std::cerr << "OscCurve's fHist is associated with a directory. How did that happen?" << std::endl;
      }
    }

    HistCache::Delete(fHist);
  }

  //----------------------------------------------------------------------
  OscCurve::OscCurve(const OscCurve& rhs)
  {
    DontAddDirectory guard;

    assert(rhs.fHist);
    fHist = HistCache::Copy(rhs.fHist);
  }

  //----------------------------------------------------------------------
  OscCurve& OscCurve::operator=(const OscCurve& rhs)
  {
    if(&rhs == this) return *this;

    DontAddDirectory guard;

    HistCache::Delete(fHist);
    assert(rhs.fHist);
    fHist = HistCache::Copy(rhs.fHist);
    return *this;
  }

  //----------------------------------------------------------------------
  TH1D* OscCurve::ToTH1(bool title) const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    TH1D* ret = HistCache::Copy(fHist);

    if(title){
      // Don't do this work unless it's explicitly requested
      std::map<int, std::string> nus;
      nus[12] = nus[-12] = "e";
      nus[14] = nus[-14] = "#mu";
      nus[16] = nus[-16] = "#tau";
      nus[0] = "active";
      const std::string nu = (fFrom > 0) ? "#nu" : "#bar#nu";

      ret->SetTitle((nu+"_{"+nus[fFrom]+"}#rightarrow"+nu+"_{"+nus[fTo]+"}").c_str());
    }

    return ret;
  }
}
