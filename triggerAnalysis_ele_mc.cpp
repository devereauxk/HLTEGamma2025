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
#include <TAxis.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "pinchun/include/CommandLine.h"
#include "triggerAnalysis_ele_mc.C"

namespace {

const std::string kDefaultZToEEInput =
  "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/"
  "JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/"
  "crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/";

void PrintUsage(const char *self)
{
  std::cout
    << "Usage: " << self << " [options]\n"
    << "  --inputForest <path>   Forest directory\n"
    << "  --inputHLT <path>      Optional HLT directory (default: empty)\n"
    << "  --outputBase <name>    Output tag for plots and ROOT file\n"
    << "  --nametag <label>      Plot label\n"
    << "  --matchToSelf <bool>   Match reco electrons to HLT objects in same forest\n"
    << "  --nfiles <int>         Number of HiForestMiniAOD_<N>.root files to read (-1 for glob)\n"
    << "  --minHiBin <float>     Minimum hiBin\n"
    << "  --maxHiBin <float>     Maximum hiBin\n";
}

}  // namespace

int main(int argc, char *argv[])
{
  CommandLine cl(argc, argv);

  if (cl.GetBool("help", false) || cl.GetBool("h", false)) {
    PrintUsage(argv[0]);
    return 0;
  }

  const std::string inputForest = cl.Get("inputForest", kDefaultZToEEInput);
  const std::string inputHLT = cl.Get("inputHLT", "");
  const std::string outputBase = cl.Get("outputBase", "MCZee0_100_2023PbPbcuts");
  const std::string nametag = cl.Get("nametag", "2025 Pythia8+Hydjet Z->EE");
  const bool matchToSelf = cl.GetBool("matchToSelf", true);
  const int nfiles = cl.GetInt("nfiles", 10);
  const float minHiBin = cl.GetDouble("minHiBin", 0.0);
  const float maxHiBin = cl.GetDouble("maxHiBin", 200.0);

  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("output", kTRUE);

  triggerAnalysis_ele_mc(
    inputForest,
    inputHLT,
    outputBase,
    nametag,
    matchToSelf,
    nfiles,
    minHiBin,
    maxHiBin
  );

  return 0;
}
