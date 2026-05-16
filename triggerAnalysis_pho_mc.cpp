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
#include "triggerAnalysis_pho_mc.C"

int main(int argc, char *argv[])
{
  CommandLine cl(argc, argv);

  const std::string inputForest = cl.Get("inputForest",
    "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2026MC/"
    "PhotonQCD_pTHatMin20_HydjetEmbedded_1610pre3/"
    "crab_Run3_PbPb_2026MC_ZToEE/260514_202916/0000/");
  const std::string inputEmulation  = cl.Get("inputEmulation", "");
  const std::string outputBase      = cl.Get("outputBase", "MCQCDPhoton0_100");
  const std::string nametag         = cl.Get("nametag", "2026 Pythia8+Hydjet QCD Photon");
  const bool        matchToForest   = cl.GetBool("matchToForest", true);
  const bool        matchToEmulation = cl.GetBool("matchToEmulation", false);
  const int         nfiles          = cl.GetInt("nfiles", -1);
  const float       minHiBin        = cl.GetDouble("minHiBin", 0.0);
  const float       maxHiBin        = cl.GetDouble("maxHiBin", 200.0);
  const int         year            = cl.GetInt("year", 2026);
  const std::string plotDir         = cl.Get("plotDir", "plots");

  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(plotDir.c_str(), kTRUE);
  gSystem->mkdir("output", kTRUE);

  triggerAnalysis_pho_mc(
    inputForest, inputEmulation, outputBase, nametag,
    matchToForest, matchToEmulation,
    nfiles, minHiBin, maxHiBin, year, plotDir);

  return 0;
}
