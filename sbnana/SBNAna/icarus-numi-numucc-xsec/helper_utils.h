// Utilities...

// ROOT
#include "TPaveText.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1D.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"

// C++
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cassert>

using namespace ana;

////////////////////////////////////////////////////////////////
// Text additions
////////////////////////////////////////////////////////////////
void PrintICARUSSim()
{
  TPaveText *tpt = new TPaveText(0.63,0.91,0.9,0.96,"NB NDC");
  TText *txt = tpt->AddText("ICARUS Simulation");
  txt->SetTextColor(kRed+1);

  tpt->Draw();
}

void PrintICARUSData()
{
  TPaveText *tpt = new TPaveText(0.63,0.91,0.9,0.96,"NB NDC");
  TText *txt = tpt->AddText("ICARUS Data");
  txt->SetTextColor(kBlue+1);

  tpt->Draw();
}

void PrintPreliminary(int vertLoc=1, int horizLoc=3)
{
  double xLo(0.58), xHi(0.85), yLo(0.78), yHi(0.86);
  if ( vertLoc==0 ) { yLo=0.14; yHi=.22; }
  if ( horizLoc==0 ) { xLo=0.15; xHi=0.42; }
  else if ( horizLoc==1 ) { xLo=0.36; xHi=0.63; }
  if ( vertLoc==-1 && horizLoc==-1 ) { xLo=0.25; yLo=0.91; xHi=0.52; yHi=0.96; }
  if ( vertLoc==-2 && horizLoc==-2 ) { xLo=0.53; yLo=0.40; xHi=0.82; yHi=0.48; }
  if ( vertLoc==-3 && horizLoc==-3 ) { xLo=0.11; yLo=0.50; xHi=0.35; yHi=0.58; }

  TPaveText *tpt = new TPaveText(xLo,yLo,xHi,yHi,"NB NDC");
  tpt->SetFillStyle(0);
  tpt->SetBorderSize(0);
  tpt->SetLineColor(0);

  TText *txt = tpt->AddText("Preliminary");
  txt->SetTextColor(kGray+1);

  tpt->Draw();
}

void PrintWIP(int vertLoc=1, int horizLoc=3)
{
  double xLo(0.55), xHi(0.85), yLo(0.78), yHi(0.86);
  if ( vertLoc==0 ) { yLo=0.14; yHi=.25; }
  if ( horizLoc==0 ) { xLo=0.15; xHi=0.42; }
  else if ( horizLoc==1 ) { xLo=0.35; xHi=0.65; }
  if ( vertLoc==-1 && horizLoc==-1 ) { xLo=0.25; yLo=0.91; xHi=0.52; yHi=0.96; }
  if ( vertLoc==-2 && horizLoc==-2 ) { xLo=0.53; yLo=0.40; xHi=0.82; yHi=0.48; }
  if ( vertLoc==-3 && horizLoc==-3 ) { xLo=0.11; yLo=0.50; xHi=0.37; yHi=0.60; }

  TPaveText *tpt = new TPaveText(xLo,yLo,xHi,yHi,"NB NDC");
  tpt->SetFillStyle(0);
  tpt->SetBorderSize(0);
  tpt->SetLineColor(0);

  TText *txt = tpt->AddText("Work in Progress");
  txt->SetTextColor(kGray+1);

  tpt->Draw();
}

void PrintTypeText(int Loc=0)
{
  double xLo(0.58), xHi(0.85), yLo(0.60), yHi(0.74);
  if ( Loc==1 ) { xLo=0.14; xHi=0.35; }

  TPaveText *tpt = new TPaveText(xLo,yLo,xHi,yHi,"NB NDC");
  tpt->SetFillStyle(0);
  tpt->SetBorderSize(0);
  tpt->SetLineColor(0);

  TText *txt = tpt->AddText("All signal");
  txt->SetTextColor(kBlack);
  TText *txt2 = tpt->AddText("Selected");
  txt2->SetTextColor(kRed);

  tpt->Draw();
}

void PrintTypeTextMore(int Loc=0)
{
  double xLo(0.58), xHi(0.85), yLo(0.60), yHi(0.85);
  if ( Loc==1 ) { xLo=0.14; xHi=0.35; }

  TPaveText *tpt = new TPaveText(xLo,yLo,xHi,yHi,"NB NDC");
  tpt->SetFillStyle(0);
  tpt->SetBorderSize(0);
  tpt->SetLineColor(0);

  TText *txt = tpt->AddText("All signal");
  txt->SetTextColor(kBlack);
  TText *txt2 = tpt->AddText("UpToMuID");
  txt2->SetTextColor(kBlue);
  TText *txt3 = tpt->AddText("UpToPrID");
  txt3->SetTextColor(kGreen+2);
  TText *txt4 = tpt->AddText("Selected");
  txt4->SetTextColor(kRed);

  tpt->Draw();
}

//////////////////////
// Plot help:
//////////////////////

std::vector<TH1*> makeHistos ( const std::vector<Spectrum*> spectra, const std::vector<int> colors, const float exposure, const bool fill=false )
{
  std::vector<TH1*> histos;
  for ( unsigned int idx = 0; idx < spectra.size(); ++idx ) {
    if( idx==0 ) histos.push_back( spectra.at(idx)->ToTH1(exposure, kBlack) );
    else if( idx!=0 ) histos.push_back( spectra.at(idx)->ToTH1(exposure, colors.at(idx-1)) );
    histos.at(idx)->GetYaxis()->SetTitle( TString::Format( "Slices / %.2fe20 POT", exposure/1.0e20) );
    if( idx!=0 && fill ) histos.at(idx)->SetFillColor(colors.at(idx-1));
  }

  return histos;
}

THStack* makeStack( const std::vector<TH1*> histos, const float exposure, std::string stackName )
{
  assert(histos.size() >= 2);

  THStack *hstack = new THStack(stackName.c_str(),"");

  histos.at(1)->GetYaxis()->SetTitle( TString::Format( "Slices / %.2fe20 POT", exposure/1.0e20) );
  hstack->Add(histos.at(1));

  for ( unsigned int idx = 2; idx < histos.size(); ++idx ) {
    hstack->Add(histos.at(idx));
  }

  return hstack;
}

std::vector<float> makeIntegrals( const std::vector<Spectrum*> spectra, const float exposure )
{
  std::vector<float> integrals;

  for ( unsigned int idx = 0; idx < spectra.size(); ++idx ) {
    integrals.push_back( spectra.at(idx)->Integral(exposure) );
  }

  return integrals;
} 

void fillLegend( TLegend* theLegend, const std::vector<TH1*> histos, const std::vector<std::string> titles,
                 const std::vector<float> integrals={} )
{
    if( integrals.size()!=0 ) theLegend->SetHeader( TString::Format("Total: %.1fx10^{3}", integrals[0]/1000.) );
    for ( unsigned int idx=1; idx < histos.size(); ++idx ) {
      if ( integrals.size()!=0 )
        theLegend->AddEntry( histos.at(idx), TString::Format("%s: %.2f%%", titles.at(idx).c_str(), 100.*integrals.at(idx)/integrals.at(0)), "f");
      else 
        theLegend->AddEntry( histos.at(idx), titles.at(idx).c_str(), "l" );
    }
}

//////////////////////
// Colors:
//////////////////////

// We'll use it as Signal, Other NuMu CC, NC, Cosmics ('out-of-time'), In-time cosmic
const int colorwheel[5] = {kBlue, kGreen+2, kRed, kMagenta, kViolet-1};
const std::vector<int> colorwheelVec = {kBlue, kGreen+2, kRed, kMagenta, kViolet-1};
const std::vector<int> colorwheelVec2 = {kBlue, kCyan+1, kGreen+2, kRed, kMagenta, kViolet-1};

// order: QE, MEC, RES, DIS, COH, NonSigCC
const int colorwheel_mode[6] = {kBlue, kCyan+1, kGreen+2, kOrange+1, kGreen, kGray+2};
const std::vector<int> colorwheelVec_mode = {kBlue, kCyan+1, kGreen+2, kOrange+1, kGreen, kGray+2};

// order: true signal, numu mu wrong, numu p wrong, numu both wrong, numu non fid, numu other, nu other, cosmic, cosmic in-time
const int colorwheel_class[9] = {kBlue, kAzure+7, kRed, kGray+2, kGreen+2, kGreen, kOrange+2, kMagenta, kViolet-1};
const std::vector<int> colorwheelVec_class = {kBlue, kAzure+7, kRed, kGray+2, kGreen+2, kGreen, kOrange+2, kMagenta, kViolet-1};
