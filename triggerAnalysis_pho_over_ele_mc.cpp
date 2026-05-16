#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "include/CommandLine.h"
#include "triggerAnalysis_pho_over_ele_mc.C"

int main(int argc, char *argv[])
{
  CommandLine cl(argc, argv);

  const std::string inputForest = cl.Get("inputForest",
    "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/"
    "JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/"
    "crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/");
  const std::string outputBase  = cl.Get("outputBase", "MCphoOverEle0_100");
  const std::string nametag     = cl.Get("nametag", "2025 Pythia8+Hydjet Z->EE");
  const int         nfiles      = cl.GetInt("nfiles", -1);
  const float       minHiBin    = cl.GetDouble("minHiBin", 0.0);
  const float       maxHiBin    = cl.GetDouble("maxHiBin", 200.0);
  const int         year        = cl.GetInt("year", 2025);
  const std::string plotDir     = cl.Get("plotDir", "plots");

  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(plotDir.c_str(), kTRUE);
  gSystem->mkdir("output", kTRUE);

  triggerAnalysis_pho_over_ele_mc(
    inputForest, outputBase, nametag,
    nfiles, minHiBin, maxHiBin, year, plotDir);

  return 0;
}
