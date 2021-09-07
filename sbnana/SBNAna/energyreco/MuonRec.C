// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

// ---- VARS -----
//const Var kErecMuon([](const caf::SRSliceProxy* slc) -> float {
const Var kErecMuon([](const caf::SRSliceProxy* slc) -> int {
  float maxlen = 0;
  int index = 0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;
  }
  for (auto const& trk : slc->reco.trk) {
    i += 1;
    if (trk.len == maxlen)
      index = i;
  }
  return index;
});

  //for (auto i=0; i<slc->reco.trk.size(); i++) {
    //if (slc->reco.trk[i].len > maxlen)
      //maxlen = slc->reco.trk[i].len;}

  //for (auto const& trk : slc->reco.trk) { 
    //if (trk.len == len)
      //index = slc;
  //return slc;

// ---- CUTS -----

// one file
void MuonRec(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/caf-01608068-6435-4a36-93b5-29ead574d963.root")

// twenty-one files
// //void Muon(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/*.root")

// all files
// void Muon(const std::string inputName = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus")

{
  SpectrumLoader loader(inputName);

  // ---- SPECTRA -----
  // A spectrum is a histogram with associated POT information
  const Binning bins1 = Binning::Simple(40, 0, 200);

  // HistAxis(Title, Binning, Var)
  const HistAxis axErecMuon("Reconstructed Muon Energy ()", bins1, kErecMuon);

  // Spectrum(Spectrumloader, HistAxis, Cut)
  Spectrum sErecMuon(loader, axErecMuon, kNoCut);

  // This is the call that actually fills in the spectrum
  loader.Go();

  // ---- DRAW -----
  // For plotting purposes we can convert spectra to a TH1
  TCanvas* c1 = new TCanvas("c1", "c1");
  TH1* hErecMuon = sErecMuon.ToTH1(6.6e20);
  hErecMuon->SetLineColor(kGreen);
  hErecMuon->Draw("hist");

  //Drawing a legend because we are not hobos
  TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg->AddEntry(hErecMuon, "Reconstructed Muon Energy ()", "l");
  leg->Draw();
  c1->Print("EMuon.png");
}
