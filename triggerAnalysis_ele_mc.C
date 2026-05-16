#include <algorithm>
#include <array>
#include <dirent.h>

using namespace std;

namespace {

constexpr size_t kNL1Paths = 3;
constexpr size_t kNHLTPaths = 5;

vector<string> collectRootFiles(const string &directory)
{
  vector<string> files;

  DIR *dir = opendir(directory.c_str());
  if (!dir) return files;

  while (dirent *entry = readdir(dir)) {
    string fileName = entry->d_name;
    if (fileName.size() >= 5 && fileName.substr(fileName.size() - 5) == ".root") {
      files.push_back(directory + "/" + fileName);
    }
  }
  closedir(dir);

  sort(files.begin(), files.end());
  return files;
}

}  // namespace

Float_t PtMin = 2.;
Float_t PtMax = 110;
Int_t NPtBins = 35;

Float_t etaMin = -5.0;
Float_t etaMax = 5.0;
Float_t NEtaBins = 100;

Float_t phiMin = -TMath::Pi();
Float_t phiMax = TMath::Pi();
Float_t NPhiBins = 100;

TH1D *denom = new TH1D("denom","denom",NPtBins,PtMin,PtMax);
TH1D *denom_endcap = new TH1D("denom_endcap","denom_endcap",NPtBins,PtMin,PtMax);

TH1D *l1_7 = new TH1D("l1_7","l1_7",NPtBins,PtMin,PtMax);
TH1D *l1_15 = new TH1D("l1_15","l1_15",NPtBins,PtMin,PtMax);
TH1D *l1_21 = new TH1D("l1_21","l1_21",NPtBins,PtMin,PtMax);

TH1D *l1_7_endcap = new TH1D("l1_7_endcap","l1_7_endcap",NPtBins,PtMin,PtMax);
TH1D *l1_15_endcap = new TH1D("l1_15_endcap","l1_15_endcap",NPtBins,PtMin,PtMax);
TH1D *l1_21_endcap = new TH1D("l1_21_endcap","l1_21_endcap",NPtBins,PtMin,PtMax);

TH1D *rl1_7 = new TH1D("rl1_7","rl1_7",NPtBins,PtMin,PtMax);
TH1D *rl1_15 = new TH1D("rl1_15","rl1_15",NPtBins,PtMin,PtMax);
TH1D *rl1_21 = new TH1D("rl1_21","rl1_21",NPtBins,PtMin,PtMax);

TH1D *rl1_7_endcap = new TH1D("rl1_7_endcap","rl1_7_endcap",NPtBins,PtMin,PtMax);
TH1D *rl1_15_endcap = new TH1D("rl1_15_endcap","rl1_15_endcap",NPtBins,PtMin,PtMax);
TH1D *rl1_21_endcap = new TH1D("rl1_21_endcap","rl1_21_endcap",NPtBins,PtMin,PtMax);

TH1D *num_15 = new TH1D("num_15","num_15",NPtBins,PtMin,PtMax);
TH1D *num_15_endcap = new TH1D("num_15_endcap","num_15_endcap",NPtBins,PtMin,PtMax);
TH1D *r_15 = new TH1D("r_15","r_15",NPtBins,PtMin,PtMax);
TH1D *r_15_endcap = new TH1D("r_15_endcap","r_15_endcap",NPtBins,PtMin,PtMax);

TH1D *num_20 = new TH1D("num_20","num_20",NPtBins,PtMin,PtMax);
TH1D *num_20_endcap = new TH1D("num_20_endcap","num_20_endcap",NPtBins,PtMin,PtMax);
TH1D *r_20 = new TH1D("r_20","r_20",NPtBins,PtMin,PtMax);
TH1D *r_20_endcap = new TH1D("r_20_endcap","r_20_endcap",NPtBins,PtMin,PtMax);

TH1D *num_30 = new TH1D("num_30","num_30",NPtBins,PtMin,PtMax);
TH1D *num_30_endcap = new TH1D("num_30_endcap","num_30_endcap",NPtBins,PtMin,PtMax);
TH1D *r_30 = new TH1D("r_30","r_30",NPtBins,PtMin,PtMax);
TH1D *r_30_endcap = new TH1D("r_30_endcap","r_30_endcap",NPtBins,PtMin,PtMax);

TH1D *num_40 = new TH1D("num_40","num_40",NPtBins,PtMin,PtMax);
TH1D *num_40_endcap = new TH1D("num_40_endcap","num_40_endcap",NPtBins,PtMin,PtMax);
TH1D *r_40 = new TH1D("r_40","r_40",NPtBins,PtMin,PtMax);
TH1D *r_40_endcap = new TH1D("r_40_endcap","r_40_endcap",NPtBins,PtMin,PtMax);

TH1D *num_50 = new TH1D("num_50","num_50",NPtBins,PtMin,PtMax);
TH1D *num_50_endcap = new TH1D("num_50_endcap","num_50_endcap",NPtBins,PtMin,PtMax);
TH1D *r_50 = new TH1D("r_50","r_50",NPtBins,PtMin,PtMax);
TH1D *r_50_endcap = new TH1D("r_50_endcap","r_50_endcap",NPtBins,PtMin,PtMax);

//TH1D *h_hiBin = new TH1D("h_hiBin","h_hiBin",250,0,250);
//TH2D *h2_missedEtaPhi = new TH2D("h2_missedEtaPhi","h2_missedEtaPhi",NPhiBins,phiMin,phiMax,NEtaBins,etaMin,etaMax);



float deltaPhi(float phi1, float phi2)
{
  float dPhi = phi1 - phi2;
  while (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return dPhi;
}


// ==========================================================
// main: trigger turn on curves for HLT electrons in MC
// ==========================================================

void triggerAnalysis_ele_mc(

    // 2025 ZtoEE sample
    string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/",
    string inputEmulation = "",
    string output_base = "MCZee0_100_2023PbPbcuts",
    string nametag = "2025 Pythia8+Hydjet Z->EE",
    bool matchToForest = true,
    bool matchToEmulation = false,

    // 2025 ZtoEE sample, embedding with thrOverEEE, thrOverEEB = 0.5 menu
    //string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/",
    //string inputHLT = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/emulation/fdamas_151X_V5_thrOverEEE_0p5/",
    //string output_base = "MCZee0_100_thrOverEEE_0p5",
    
    // 2025 JpsiToEE sample
    //string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000/",
    //string inputHLT = "",
    //string output_base = "MCJpsiToEE0_100_2023PbPbcuts",
    //string nametag = "2025 Pythia8+Hydjet J/psi->EE",

    // 2024 ZtoEE sample
    //string inputForest = "/eos/cms/store/group/phys_heavyions/prdas/EGamma/Run3_PbPb_2024_MC/Ze10e10/Embedded/Pythia8_Embedded_Ze10e10_TuneCP5_2024/crab_20241026_145435/241026_125438/0000",
    //string inputHLT = "",
    //string output_base = "2024MCZee0_100",

    // 2024 ZtoEE sample, noembedding
    //string inputForest = "/eos/cms/store/group/phys_heavyions/prdas/EGamma/Run3_PbPb_2024_MC/Ze10e10/NoEmbedding/Pythia8_Embedded_Ze10e10_TuneCP5_2024/crab_20241026_145054/241026_125057/0000",
    //string inputHLT = "",
    //string output_base = "2024MCZeeNoembedding0_100",

    // 2024 ppref, ZtoEE
    //string inputForest = "/eos/cms/store/group/phys_heavyions/prdas/EGamma/Run3_ppref_2024_MC/Ze10e10/Pythia8_ppRef_Ze10e10_TuneCP5_2024/crab_20241026_131030/241026_111033/0000",
    //string inputHLT = "",
    //string output_base = "2024MCpprefZee",

    int nfiles = 10,
    float minHiBin = 0.0,
    float maxHiBin = 200.0,
    int year = 2025,
    string plotDir = "plots"
  ){

  const string hlt_ele15 = (year == 2026) ? "HLT_HIEle15Gsf_v18" : "HLT_HIEle15Gsf_v16";
  const string hlt_ele20 = (year == 2026) ? "HLT_HIEle20Gsf_v18" : "HLT_HIEle20Gsf_v16";
  const string hlt_ele30 = (year == 2026) ? "HLT_HIEle30Gsf_v18" : "HLT_HIEle30Gsf_v16";
  const string hlt_ele40 = (year == 2026) ? "HLT_HIEle40Gsf_v18" : "HLT_HIEle40Gsf_v16";
  const string hlt_ele50 = (year == 2026) ? "HLT_HIEle50Gsf_v18" : "HLT_HIEle50Gsf_v16";

  std::cout << "running triggerAnalysis_ele_mc()" << std::endl;
  std::cout << "input forest directory = " << inputForest  << std::endl;
  std::cout << "input emulation dir    = " << inputEmulation  << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;
  std::cout << "matchToForest          = " << matchToForest << std::endl;
  std::cout << "matchToEmulation       = " << matchToEmulation << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");
  vector<string> HLTObjectTreeNames = {
    "hltobject/HLT_HIEle15Gsf_v",
    "hltobject/HLT_HIEle20Gsf_v",
    "hltobject/HLT_HIEle30Gsf_v",
    "hltobject/HLT_HIEle40Gsf_v",
    "hltobject/HLT_HIEle50Gsf_v"
  };
  
  std::cout<< "Adding input files...";

  // forest files
  if (nfiles == -1) {
    string this_forest_name = inputForest + "/HiForestMiniAOD_*";
    HltTree   ->Add(this_forest_name.c_str());
    EventTree ->Add(this_forest_name.c_str());
    HiTree    ->Add(this_forest_name.c_str());
  }
  else {
    for(int i = 1; i <= nfiles; i++) {
      string this_forest_name = inputForest + "/HiForestMiniAOD_" + to_string(i) + ".root";
      HltTree   ->Add(this_forest_name.c_str());
      EventTree ->Add(this_forest_name.c_str());
      HiTree    ->Add(this_forest_name.c_str());
    }
  }

  std::cout << "done" << std::endl;

  std::cout << "Setting Event, lumi, and run branchAdresses...";

  // ================= Hi Tree =================
  ULong64_t       evt;
  UInt_t          lumi;
  UInt_t          run;
  Int_t           hiBin;
  Float_t         hiHF;
  Float_t         pthat_weight;
  Float_t         vz;
  
  HiTree->SetBranchStatus("*",0);
  HiTree->SetBranchStatus("evt",1);
  HiTree->SetBranchStatus("lumi",1);
  HiTree->SetBranchStatus("run",1);
  HiTree->SetBranchStatus("hiBin",1);
  HiTree->SetBranchStatus("hiHF",1);
  HiTree->SetBranchStatus("weight",1);
  HiTree->SetBranchStatus("vz",1);

  HiTree->SetBranchAddress("evt", &evt);
  HiTree->SetBranchAddress("lumi", &lumi);
  HiTree->SetBranchAddress("run", &run);
  HiTree->SetBranchAddress("hiBin", &hiBin);
  HiTree->SetBranchAddress("hiHF", &hiHF);
  HiTree->SetBranchAddress("weight", &pthat_weight);
  HiTree->SetBranchAddress("vz", &vz);


  // ================= HLT Trees =================
  // forest
  Int_t           L1_SingleEG7;
  Int_t           L1_SingleEG15;
  Int_t           L1_SingleEG21;
  Int_t           HLT_HIEle15Gsf;
  Int_t           HLT_HIEle20Gsf;
  Int_t           HLT_HIEle30Gsf;
  Int_t           HLT_HIEle40Gsf;
  Int_t           HLT_HIEle50Gsf;
  ULong64_t       forest_Event;
  Int_t           forest_LumiBlock, forest_Run;

  HltTree->SetBranchStatus("*",0);     // disable all branches
  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  HltTree->SetBranchStatus("L1_SingleEG7_BptxAND", 1);
  HltTree->SetBranchStatus("L1_SingleEG15_BptxAND", 1);
  HltTree->SetBranchStatus("L1_SingleEG21_BptxAND", 1);

  HltTree->SetBranchStatus(hlt_ele15.c_str(), 1);
  HltTree->SetBranchStatus(hlt_ele20.c_str(), 1);
  HltTree->SetBranchStatus(hlt_ele30.c_str(), 1);
  HltTree->SetBranchStatus(hlt_ele40.c_str(), 1);
  HltTree->SetBranchStatus(hlt_ele50.c_str(), 1);

  HltTree->SetBranchAddress("Event", &forest_Event);
  HltTree->SetBranchAddress("LumiBlock", &forest_LumiBlock);
  HltTree->SetBranchAddress("Run", &forest_Run);

  HltTree->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7);
  HltTree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15);
  HltTree->SetBranchAddress("L1_SingleEG21_BptxAND", &L1_SingleEG21);

  HltTree->SetBranchAddress(hlt_ele15.c_str(), &HLT_HIEle15Gsf);
  HltTree->SetBranchAddress(hlt_ele20.c_str(), &HLT_HIEle20Gsf);
  HltTree->SetBranchAddress(hlt_ele30.c_str(), &HLT_HIEle30Gsf);
  HltTree->SetBranchAddress(hlt_ele40.c_str(), &HLT_HIEle40Gsf);
  HltTree->SetBranchAddress(hlt_ele50.c_str(), &HLT_HIEle50Gsf);

  vector<float>           *Ele15Gsf_pt = nullptr;
  vector<float>           *Ele15Gsf_eta = nullptr;
  vector<float>           *Ele15Gsf_phi = nullptr;
  vector<float>           *Ele20Gsf_pt = nullptr;
  vector<float>           *Ele20Gsf_eta = nullptr;
  vector<float>           *Ele20Gsf_phi = nullptr;
  vector<float>           *Ele30Gsf_pt = nullptr;
  vector<float>           *Ele30Gsf_eta = nullptr;
  vector<float>           *Ele30Gsf_phi = nullptr;
  vector<float>           *Ele40Gsf_pt = nullptr;
  vector<float>           *Ele40Gsf_eta = nullptr;
  vector<float>           *Ele40Gsf_phi = nullptr;
  vector<float>           *Ele50Gsf_pt = nullptr;
  vector<float>           *Ele50Gsf_eta = nullptr;
  vector<float>           *Ele50Gsf_phi = nullptr;

  vector<vector<float>*>   HLTObject_pt = {Ele15Gsf_pt, Ele20Gsf_pt, Ele30Gsf_pt, Ele40Gsf_pt, Ele50Gsf_pt};
  vector<vector<float>*>   HLTObject_eta = {Ele15Gsf_eta, Ele20Gsf_eta, Ele30Gsf_eta, Ele40Gsf_eta, Ele50Gsf_eta};
  vector<vector<float>*>   HLTObject_phi = {Ele15Gsf_phi, Ele20Gsf_phi, Ele30Gsf_phi, Ele40Gsf_phi, Ele50Gsf_phi};
  vector<TTree*> currentHLTObjectTrees(HLTObjectTreeNames.size(), nullptr);
  int currentHLTTreeNumber = -1;

  // ===== Emulation TChain setup (only used when matchToEmulation=true) =====
  TChain *EmulationHltTree = nullptr;
  Int_t HLT_emul_HIEle15Gsf = 0, HLT_emul_HIEle20Gsf = 0, HLT_emul_HIEle30Gsf = 0;
  Int_t HLT_emul_HIEle40Gsf = 0, HLT_emul_HIEle50Gsf = 0;

  if (matchToEmulation) {
    if (inputEmulation.empty()) {
      std::cerr << "matchToEmulation requires --inputEmulation." << std::endl;
      return;
    }
    const vector<string> emulationFiles = collectRootFiles(inputEmulation);
    if (emulationFiles.empty()) {
      std::cerr << "No emulation ROOT files found in " << inputEmulation << std::endl;
      return;
    }
    EmulationHltTree = new TChain("hltanalysis/HltTree");
    for (const string &f : emulationFiles) EmulationHltTree->Add(f.c_str());

    EmulationHltTree->SetBranchStatus("*", 0);
    EmulationHltTree->SetBranchStatus(hlt_ele15.c_str(), 1);
    EmulationHltTree->SetBranchStatus(hlt_ele20.c_str(), 1);
    EmulationHltTree->SetBranchStatus(hlt_ele30.c_str(), 1);
    EmulationHltTree->SetBranchStatus(hlt_ele40.c_str(), 1);
    EmulationHltTree->SetBranchStatus(hlt_ele50.c_str(), 1);

    EmulationHltTree->SetBranchAddress(hlt_ele15.c_str(), &HLT_emul_HIEle15Gsf);
    EmulationHltTree->SetBranchAddress(hlt_ele20.c_str(), &HLT_emul_HIEle20Gsf);
    EmulationHltTree->SetBranchAddress(hlt_ele30.c_str(), &HLT_emul_HIEle30Gsf);
    EmulationHltTree->SetBranchAddress(hlt_ele40.c_str(), &HLT_emul_HIEle40Gsf);
    EmulationHltTree->SetBranchAddress(hlt_ele50.c_str(), &HLT_emul_HIEle50Gsf);

    std::cout << "Emulation chain: " << EmulationHltTree->GetEntries() << " entries from "
              << emulationFiles.size() << " files" << std::endl;
    if (EmulationHltTree->GetEntries() != HltTree->GetEntries()) {
      std::cerr << "Warning: emulation (" << EmulationHltTree->GetEntries()
                << ") and forest (" << HltTree->GetEntries()
                << ") entry counts differ — check alignment." << std::endl;
    }
  }


  // ================= Event Tree =================
  vector<float>   *elePt = nullptr;
  vector<float>   *eleEta = nullptr;
  vector<float>   *elePhi = nullptr;
  vector<int>     *eleCharge = nullptr;
  Int_t           nEle;
  vector<int>     *ele_genMatchedIndex = nullptr;
  vector<float>   *eleHoverE = nullptr;
  vector<float>   *eleSigmaIEtaIEta_2012 = nullptr;
  vector<float>   *eledEtaSeedAtVtx = nullptr;
  vector<float>   *eledPhiAtVtx = nullptr;
  vector<float>   *eleEoverPInv = nullptr;
  vector<int>     *eleMissHits = nullptr;
  vector<float>   *eleD0 = nullptr;
  vector<float>   *eleSCEta = nullptr;
  vector<float>   *eleSCRawEn = nullptr;
  vector<int>     *mcStatus = nullptr;
  vector<int>     *mcPID = nullptr;
  vector<int>     *mcMomPID = nullptr;
  vector<float>   *mcPt = nullptr;
  vector<float>   *mcEta = nullptr;
  vector<float>   *mcPhi = nullptr;
  vector<float>   *eleDz = nullptr;

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

  // mc branches
  EventTree->SetBranchStatus("mcStatus",1);
  EventTree->SetBranchStatus("mcPID",1);
  EventTree->SetBranchStatus("mcMomPID",1);
  EventTree->SetBranchStatus("mcPt",1);
  EventTree->SetBranchStatus("mcEta",1);
  EventTree->SetBranchStatus("mcPhi",1);

  EventTree->SetBranchStatus("nEle",1);
  EventTree->SetBranchStatus("elePt",1);
  EventTree->SetBranchStatus("eleEta",1);
  EventTree->SetBranchStatus("eleHoverE",1);
  EventTree->SetBranchStatus("elePhi",1);
  EventTree->SetBranchStatus("eleCharge",1);
  //EventTree->SetBranchStatus("ele_genMatchedIndex",1);

  EventTree->SetBranchStatus("eleSigmaIEtaIEta_2012",1);
  EventTree->SetBranchStatus("eledEtaSeedAtVtx",1);
  EventTree->SetBranchStatus("eledPhiAtVtx",1);
  EventTree->SetBranchStatus("eleEoverPInv",1);
  EventTree->SetBranchStatus("eleMissHits",1);
  EventTree->SetBranchStatus("eleD0",1);
  EventTree->SetBranchStatus("eleSCEta",1);
  EventTree->SetBranchStatus("eleSCRawEn",1);
  EventTree->SetBranchStatus("eleDz",1);

  EventTree->SetBranchAddress("mcStatus", &mcStatus);
  EventTree->SetBranchAddress("mcPID", &mcPID);
  EventTree->SetBranchAddress("mcMomPID", &mcMomPID);
  EventTree->SetBranchAddress("mcPt", &mcPt);
  EventTree->SetBranchAddress("mcEta", &mcEta);
  EventTree->SetBranchAddress("mcPhi", &mcPhi);
  //EventTree->SetBranchAddress("ele_genMatchedIndex", &ele_genMatchedIndex);

  EventTree->SetBranchAddress("nEle", &nEle);
  EventTree->SetBranchAddress("elePt", &elePt);
  EventTree->SetBranchAddress("eleEta", &eleEta);
  EventTree->SetBranchAddress("eleHoverE", &eleHoverE);
  EventTree->SetBranchAddress("elePhi", &elePhi);
  EventTree->SetBranchAddress("eleCharge",&eleCharge);

  EventTree->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
  EventTree->SetBranchAddress("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
  EventTree->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
  EventTree->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
  EventTree->SetBranchAddress("eleMissHits", &eleMissHits);
  EventTree->SetBranchAddress("eleD0", &eleD0);
  EventTree->SetBranchAddress("eleSCEta", &eleSCEta);
  EventTree->SetBranchAddress("eleSCRawEn", &eleSCRawEn);
  EventTree->SetBranchAddress("eleDz", &eleDz);
  

  std::cout << "done" << std::endl;

  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop  =================

  // loop through reco objects
  int neles = 0;
  int npass_hltmatch = 0;
  int npass_probe = 0;
  int npass_l1 = 0;
  int npass_trigger = 0;
  for (ULong64_t i_event = 0; i_event < HltTree->GetEntries(); ++i_event){

    if(i_event%(HltTree->GetEntries()/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    Long64_t localEntry = HltTree->LoadTree(i_event);
    HiTree->LoadTree(i_event);
    EventTree->LoadTree(i_event);

    HiTree->GetEntry(i_event);
    HltTree->GetEntry(i_event);
    EventTree->GetEntry(i_event);

    //h_hiBin->Fill(hiBin);

    // event cuts
    if(fabs(vz)>15.0) continue;
    if(hiBin < minHiBin || hiBin >= maxHiBin) continue;

    neles += nEle;

    // ===== load emulation pass bits by entry index (assumes 1-1 alignment with forest) =====
    Long64_t localEmulationEntry = -1;
    if (matchToEmulation) {
      localEmulationEntry = EmulationHltTree->LoadTree(i_event);
      EmulationHltTree->GetEntry(i_event);
    }

    // ===== load HLT object pt/eta/phi from forest or emulation file =====
    if (matchToForest || matchToEmulation) {
      const int currentTreeNum = matchToForest ? HltTree->GetTreeNumber()
                                               : EmulationHltTree->GetTreeNumber();
      if (currentTreeNum != currentHLTTreeNumber) {
        currentHLTTreeNumber = currentTreeNum;
        TFile *currentFile = matchToForest ? HltTree->GetCurrentFile()
                                           : EmulationHltTree->GetCurrentFile();
        if (!currentFile) {
          std::cerr << "Could not access current file for HLT object matching." << std::endl;
          return;
        }

        for (int i_hlt = 0; i_hlt < (int)HLTObjectTreeNames.size(); i_hlt++) {
          currentHLTObjectTrees[i_hlt] = dynamic_cast<TTree*>(currentFile->Get(HLTObjectTreeNames.at(i_hlt).c_str()));
          if (!currentHLTObjectTrees[i_hlt]) {
            std::cerr << "Missing HLT object tree " << HLTObjectTreeNames.at(i_hlt)
                      << " in file " << currentFile->GetName() << std::endl;
            return;
          }

          currentHLTObjectTrees[i_hlt]->SetBranchStatus("*", 0);
          currentHLTObjectTrees[i_hlt]->SetBranchStatus("pt", 1);
          currentHLTObjectTrees[i_hlt]->SetBranchStatus("eta", 1);
          currentHLTObjectTrees[i_hlt]->SetBranchStatus("phi", 1);
          currentHLTObjectTrees[i_hlt]->SetBranchAddress("pt", &HLTObject_pt.at(i_hlt));
          currentHLTObjectTrees[i_hlt]->SetBranchAddress("eta", &HLTObject_eta.at(i_hlt));
          currentHLTObjectTrees[i_hlt]->SetBranchAddress("phi", &HLTObject_phi.at(i_hlt));
        }
      }

      const Long64_t localHLTEntry = matchToForest ? localEntry : localEmulationEntry;
      for (auto tree : currentHLTObjectTrees) tree->GetEntry(localHLTEntry);
    }

    // ====== single pass over electrons: find leading quality object per region
    //        and the best-matched quality object per HLT path per region ======
    //
    //  Denominator : leading quality electron pT, filled when L1 fires (no HLT requirement)
    //  Numerator   : pT of the quality electron with smallest dR to that path's HLT objects,
    //                filled when the HLT path fires and a match exists (dR < 0.1)

    int   i_lead_barrel = -1;  float lead_barrel_pt = -1.f;
    int   i_lead_endcap = -1;  float lead_endcap_pt = -1.f;

    for(Int_t i_track = 0; i_track < (Int_t)elePt->size(); i_track++){

      bool isBarrel = false;
      bool isEndcap = false;

      // WP95 cut-based ID (max threshold across centrality bins)
      if (abs(eleEta->at(i_track)) < 1.442) {
        if (eleMissHits->at(i_track) > 2)                 continue;
        if (abs(eleD0->at(i_track)) > 0.05)               continue;
        if (abs(eleDz->at(i_track)) > 0.10)               continue;
        if (eleSigmaIEtaIEta_2012->at(i_track) > 0.0131)  continue;
        if (abs(eledEtaSeedAtVtx->at(i_track)) > 0.00457) continue;
        if (abs(eledPhiAtVtx->at(i_track)) > 0.0963)      continue;
        if (eleEoverPInv->at(i_track) > 0.421)            continue;
        if (eleHoverE->at(i_track) > 0.156)               continue;
        isBarrel = true;
      }
      else if (abs(eleEta->at(i_track)) > 1.556 && abs(eleEta->at(i_track)) < 2.1) {
        if (eleMissHits->at(i_track) > 3)                 continue;
        if (abs(eleD0->at(i_track)) > 0.10)               continue;
        if (abs(eleDz->at(i_track)) > 0.20)               continue;
        if (eleSigmaIEtaIEta_2012->at(i_track) > 0.0382)  continue;
        if (abs(eledEtaSeedAtVtx->at(i_track)) > 0.00881) continue;
        if (abs(eledPhiAtVtx->at(i_track)) > 0.264)       continue;
        if (eleEoverPInv->at(i_track) > 0.146)            continue;
        if (eleHoverE->at(i_track) > 0.178)               continue;
        isEndcap = true;
      }
      else continue;

      if (isBarrel && elePt->at(i_track) > lead_barrel_pt) {
        lead_barrel_pt = elePt->at(i_track);
        i_lead_barrel  = i_track;
      }
      if (isEndcap && elePt->at(i_track) > lead_endcap_pt) {
        lead_endcap_pt = elePt->at(i_track);
        i_lead_endcap  = i_track;
      }
    }

    // skip events with no quality electron in either region
    if (i_lead_barrel == -1 && i_lead_endcap == -1) continue;
    npass_probe++;

    // check whether the leading electron in each region has a dR match to each HLT path
    array<bool, kNHLTPaths> lead_barrel_matched; lead_barrel_matched.fill(false);
    array<bool, kNHLTPaths> lead_endcap_matched; lead_endcap_matched.fill(false);
    if (matchToForest || matchToEmulation) {
      for (int i_hlt = 0; i_hlt < (int)currentHLTObjectTrees.size(); i_hlt++) {
        if (i_lead_barrel != -1) {
          for (int j = 0; j < (int)HLTObject_pt.at(i_hlt)->size(); j++) {
            float dEta = eleEta->at(i_lead_barrel) - HLTObject_eta.at(i_hlt)->at(j);
            float dPhi = deltaPhi(elePhi->at(i_lead_barrel), HLTObject_phi.at(i_hlt)->at(j));
            if (sqrt(dEta*dEta + dPhi*dPhi) < 0.1) { lead_barrel_matched[i_hlt] = true; break; }
          }
        }
        if (i_lead_endcap != -1) {
          for (int j = 0; j < (int)HLTObject_pt.at(i_hlt)->size(); j++) {
            float dEta = eleEta->at(i_lead_endcap) - HLTObject_eta.at(i_hlt)->at(j);
            float dPhi = deltaPhi(elePhi->at(i_lead_endcap), HLTObject_phi.at(i_hlt)->at(j));
            if (sqrt(dEta*dEta + dPhi*dPhi) < 0.1) { lead_endcap_matched[i_hlt] = true; break; }
          }
        }
      }
    }

    array<int, kNL1Paths>  L1PassMask  = {L1_SingleEG7, L1_SingleEG15, L1_SingleEG21};
    array<int, kNHLTPaths> HLTPassMask = {HLT_HIEle15Gsf, HLT_HIEle20Gsf, HLT_HIEle30Gsf, HLT_HIEle40Gsf, HLT_HIEle50Gsf};
    if (matchToEmulation) {
      HLTPassMask = {HLT_emul_HIEle15Gsf, HLT_emul_HIEle20Gsf, HLT_emul_HIEle30Gsf, HLT_emul_HIEle40Gsf, HLT_emul_HIEle50Gsf};
    }

    float weight = pthat_weight;

    // ====== fill barrel histograms (at most once per event) ======
    if (i_lead_barrel != -1) {

      float denomPt = elePt->at(i_lead_barrel);
      denom->Fill(denomPt, weight);

      if(L1PassMask[0]) {
        l1_7->Fill(denomPt, weight);
        npass_l1++;
        if(HLTPassMask[0] && lead_barrel_matched[0])
          num_15->Fill(denomPt, weight);
      }
      if(L1PassMask[1]) {
        l1_15->Fill(denomPt, weight);
        if(HLTPassMask[1] && lead_barrel_matched[1]) {
          num_20->Fill(denomPt, weight);
          npass_trigger++;
        }
        if(HLTPassMask[2] && lead_barrel_matched[2])
          num_30->Fill(denomPt, weight);
      }
      if(L1PassMask[2]) {
        l1_21->Fill(denomPt, weight);
        if(HLTPassMask[3] && lead_barrel_matched[3])
          num_40->Fill(denomPt, weight);
        if(HLTPassMask[4] && lead_barrel_matched[4])
          num_50->Fill(denomPt, weight);
      }
    }

    // ====== fill endcap histograms (at most once per event) ======
    if (i_lead_endcap != -1) {

      float denomPt = elePt->at(i_lead_endcap);
      denom_endcap->Fill(denomPt, weight);

      if(L1PassMask[0]) {
        l1_7_endcap->Fill(denomPt, weight);
        if(HLTPassMask[0] && lead_endcap_matched[0])
          num_15_endcap->Fill(denomPt, weight);
      }
      if(L1PassMask[1]) {
        l1_15_endcap->Fill(denomPt, weight);
        if(HLTPassMask[1] && lead_endcap_matched[1])
          num_20_endcap->Fill(denomPt, weight);
        if(HLTPassMask[2] && lead_endcap_matched[2])
          num_30_endcap->Fill(denomPt, weight);
      }
      if(L1PassMask[2]) {
        l1_21_endcap->Fill(denomPt, weight);
        if(HLTPassMask[3] && lead_endcap_matched[3])
          num_40_endcap->Fill(denomPt, weight);
        if(HLTPassMask[4] && lead_endcap_matched[4])
          num_50_endcap->Fill(denomPt, weight);
      }
    }

  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIEle20Gsf results ... " << std::endl;
  std::cout << "neles = " << neles << std::endl;
  std::cout << "npass_hltmatch = " << npass_hltmatch  << std::endl;
  std::cout << "npass_probe = " << npass_probe  << std::endl;
  std::cout << "npass_l1 = " << npass_l1  << std::endl;
  std::cout << "npass_trigger = " << npass_trigger  << std::endl;

  /*
  r_20->Divide(num_20,denom,1,1,"B");
  r_30->Divide(num_30,denom,1,1,"B");
  r_40->Divide(num_40,denom,1,1,"B");
  r_50->Divide(num_50,denom,1,1,"B");
  */

  r_15->Divide(num_15,l1_7,1,1,"B");
  r_20->Divide(num_20,l1_15,1,1,"B");
  r_30->Divide(num_30,l1_15,1,1,"B");
  r_40->Divide(num_40,l1_21,1,1,"B");
  r_50->Divide(num_50,l1_21,1,1,"B"); 

  rl1_7->Divide(l1_7,denom,1,1,"B");
  rl1_15->Divide(l1_15,denom,1,1,"B");
  rl1_21->Divide(l1_21,denom,1,1,"B");
  
  r_15_endcap->Divide(num_15_endcap,l1_7_endcap,1,1,"B");
  r_20_endcap->Divide(num_20_endcap,l1_15_endcap,1,1,"B");
  r_30_endcap->Divide(num_30_endcap,l1_15_endcap,1,1,"B");
  r_40_endcap->Divide(num_40_endcap,l1_21_endcap,1,1,"B");
  r_50_endcap->Divide(num_50_endcap,l1_21_endcap,1,1,"B");

  r_15->SetLineColor(kOrange+7);
  r_20->SetLineColor(kRed-4);
  r_30->SetLineColor(kBlue-4);
  r_40->SetLineColor(kGreen+2);
  r_50->SetLineColor(kMagenta-9);
  r_15_endcap->SetLineColor(kOrange+7);
  r_20_endcap->SetLineColor(kRed-4);
  r_30_endcap->SetLineColor(kBlue-4);
  r_40_endcap->SetLineColor(kGreen+2);
  r_50_endcap->SetLineColor(kMagenta-9);

  r_15->SetMarkerColor(kOrange+7);
  r_20->SetMarkerColor(kRed-4);
  r_30->SetMarkerColor(kBlue-4);
  r_40->SetMarkerColor(kGreen+2);
  r_50->SetMarkerColor(kMagenta-9);
  r_15_endcap->SetMarkerColor(kOrange+7);
  r_20_endcap->SetMarkerColor(kRed-4);
  r_30_endcap->SetMarkerColor(kBlue-4);
  r_40_endcap->SetMarkerColor(kGreen+2);
  r_50_endcap->SetMarkerColor(kMagenta-9);

  double line_width = 1.8;
  r_15->SetLineWidth(line_width);
  r_20->SetLineWidth(line_width);
  r_30->SetLineWidth(line_width);
  r_40->SetLineWidth(line_width);
  r_50->SetLineWidth(line_width);
  r_15_endcap->SetLineWidth(line_width);
  r_20_endcap->SetLineWidth(line_width);
  r_30_endcap->SetLineWidth(line_width);
  r_40_endcap->SetLineWidth(line_width);
  r_50_endcap->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_15->SetMarkerSize(marker_size);
  r_20->SetMarkerSize(marker_size);
  r_30->SetMarkerSize(marker_size);
  r_40->SetMarkerSize(marker_size);
  r_50->SetMarkerSize(marker_size);
  r_15_endcap->SetMarkerSize(marker_size);
  r_20_endcap->SetMarkerSize(marker_size);
  r_30_endcap->SetMarkerSize(marker_size);
  r_40_endcap->SetMarkerSize(marker_size);
  r_50_endcap->SetMarkerSize(marker_size);
    
  r_15->SetMarkerStyle(29);
  r_20->SetMarkerStyle(20);
  r_30->SetMarkerStyle(21);
  r_40->SetMarkerStyle(22);
  r_50->SetMarkerStyle(23);
  r_15_endcap->SetMarkerStyle(29);
  r_20_endcap->SetMarkerStyle(20);
  r_30_endcap->SetMarkerStyle(21);
  r_40_endcap->SetMarkerStyle(22);
  r_50_endcap->SetMarkerStyle(23);

  rl1_7->SetLineColor(kRed-4);
  rl1_15->SetLineColor(kBlue-4);
  rl1_21->SetLineColor(kGreen+2);

  rl1_7->SetMarkerColor(kRed-4);
  rl1_15->SetMarkerColor(kBlue-4);
  rl1_21->SetMarkerColor(kGreen+2);

  rl1_7->SetLineWidth(line_width);
  rl1_15->SetLineWidth(line_width);
  rl1_21->SetLineWidth(line_width);

  rl1_7->SetMarkerSize(marker_size);
  rl1_15->SetMarkerSize(marker_size); 
  rl1_21->SetMarkerSize(marker_size);

  rl1_7->SetMarkerStyle(20);
  rl1_15->SetMarkerStyle(21);
  rl1_21->SetMarkerStyle(22);

  // =============== Fig HLT (no L1 emulation) BARREL ==================

  TCanvas *c1 = new TCanvas("c1","c1",700,600);
  c1->cd();
  TPad *p1 = new TPad("p1","p1",0,0,1,1);
  p1->SetLeftMargin(0.13);
  p1->SetBottomMargin(0.14);
  p1->Draw();
  p1->cd();
  r_20->SetTitle("");
  r_20->SetStats(0);
  r_20->GetXaxis()->SetTitleSize(0.05);
  r_20->GetYaxis()->SetTitleSize(0.05);
  r_20->GetYaxis()->SetRangeUser(0,1.05);
  r_20->GetXaxis()->SetTitle("Reco electron #font[52]{p}_{T} [GeV]");
  r_20->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_15,"HLT_HIEle15Gsf");
  leg->AddEntry(r_20,"HLT_HIEle20Gsf");
  leg->AddEntry(r_30,"HLT_HIEle30Gsf");
  leg->AddEntry(r_40,"HLT_HIEle40Gsf");
  leg->AddEntry(r_50,"HLT_HIEle50Gsf");
  //leg->SetBorderSize(0);
  r_20->Draw();
  leg->Draw();
  r_30->Draw("same");
  r_40->Draw("same");
  r_50->Draw("same");
  r_15->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,Form("%s", nametag.c_str()));
  la->DrawLatexNDC(0.58,0.92,"Reco-HLT object matched");
  la->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("%s/TriggerEfficiency_ele_%s_barrel.png", plotDir.c_str(), output_base.c_str()));


  // =============== Fig HLT (no L1 emulation) ENCAP ==================
  TCanvas *c1_endcap = new TCanvas("c1_endcap","c1_endcap",700,600);
  c1_endcap->cd();
  TPad *p1_endcap = new TPad("p1_endcap","p1_endcap",0,0,1,1);
  p1_endcap->SetLeftMargin(0.13);
  p1_endcap->SetBottomMargin(0.14);
  p1_endcap->Draw();
  p1_endcap->cd();
  r_20_endcap->SetTitle("");
  r_20_endcap->SetStats(0);
  r_20_endcap->GetXaxis()->SetTitleSize(0.05);
  r_20_endcap->GetYaxis()->SetTitleSize(0.05);
  r_20_endcap->GetYaxis()->SetRangeUser(0,1.05);
  r_20_endcap->GetXaxis()->SetTitle("Reco electron #font[52]{p}_{T} [GeV]");
  r_20_endcap->GetYaxis()->SetTitle("Trigger efficiency");
  r_15_endcap->SetTitle("");
  r_15_endcap->SetStats(0);
  TLegend *leg_endcap = new TLegend(0.55,0.3,0.88,0.5);
  leg_endcap->AddEntry(r_15_endcap,"HLT_HIEle15Gsf");
  leg_endcap->AddEntry(r_20_endcap,"HLT_HIEle20Gsf");
  leg_endcap->AddEntry(r_30_endcap,"HLT_HIEle30Gsf");
  leg_endcap->AddEntry(r_40_endcap,"HLT_HIEle40Gsf");
  leg_endcap->AddEntry(r_50_endcap,"HLT_HIEle50Gsf");
  //leg_endcap->SetBorderSize(0);
  r_20_endcap->Draw();
  leg_endcap->Draw();
  r_30_endcap->Draw("same");
  r_40_endcap->Draw("same");
  r_50_endcap->Draw("same");
  r_15_endcap->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,Form("%s", nametag.c_str()));
  la_endcap->DrawLatexNDC(0.58,0.92,"Reco-HLT object matched");
  la_endcap->DrawLatexNDC(0.6,0.65,"Endcap: 1.556 < |#eta| < 2.1");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("%s/TriggerEfficiency_ele_%s_endcap.png", plotDir.c_str(), output_base.c_str()));


  // =============== Fig L1 ==================

  TCanvas *c2 = new TCanvas("c2","c2",700,600);
  c2->cd();
  TPad *p2 = new TPad("p2","p2",0,0,1,1);
  p2->SetLeftMargin(0.13);
  p2->SetBottomMargin(0.14);
  p2->Draw();
  p2->cd();
  rl1_7->SetTitle("");
  rl1_7->SetStats(0);
  rl1_7->GetXaxis()->SetTitleSize(0.05);
  rl1_7->GetYaxis()->SetTitleSize(0.05);
  rl1_7->GetYaxis()->SetRangeUser(0,1.05);
  rl1_7->GetXaxis()->SetTitle("Reco electron #font[52]{p}_{T} [GeV]");
  rl1_7->GetYaxis()->SetTitle("L1 Trigger efficiency");
  TLegend *leg2 = new TLegend(0.55,0.3,0.88,0.5);
  leg2->AddEntry(rl1_7,"L1_SingleEG7_BptxAND");
  leg2->AddEntry(rl1_15,"L1_SingleEG15_BptxAND");
  leg2->AddEntry(rl1_21,"L1_SingleEG21_BptxAND");
  //leg2->SetBorderSize(0);
  rl1_7->Draw();
  leg2->Draw();
  rl1_15->Draw("same");
  rl1_21->Draw("same");

  TLatex *la2 = new TLatex();
  la2->SetTextFont(42);
  la2->SetTextSize(0.03);

  la2->DrawLatexNDC(0.22,0.92,Form("%s", nametag.c_str()));
  la2->DrawLatexNDC(0.72,0.92,"L1 seed");
  la2->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la2->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c2->SaveAs(Form("%s/TriggerEfficiency_ele_%s_L1.png", plotDir.c_str(), output_base.c_str()));


  // ====== hiBin histogram — commented out ======
  //TCanvas *c3 = new TCanvas("c3","c3",700,600);
  //c3->cd();
  //TPad *p3 = new TPad("p3","p3",0,0,1,1);
  //p3->SetLeftMargin(0.13); p3->SetBottomMargin(0.14); p3->Draw(); p3->cd();
  //h_hiBin->SetTitle(""); h_hiBin->SetStats(0);
  //h_hiBin->GetXaxis()->SetTitle("hiBin"); h_hiBin->GetYaxis()->SetTitle("Counts");
  //h_hiBin->Draw();
  //c3->SaveAs(Form("%s/hiBin_ele_%s.png", plotDir.c_str(), output_base.c_str()));

  // ========== missed eta-phi plot — commented out ==========
  //TCanvas *c4 = new TCanvas("c4","c4",700,600);
  //c4->cd();
  //TPad *p4 = new TPad("p4","p4",0,0,1,1);
  //p4->SetLeftMargin(0.13); p4->SetBottomMargin(0.14); p4->Draw(); p4->cd();
  //h2_missedEtaPhi->SetTitle(""); h2_missedEtaPhi->SetStats(0);
  //h2_missedEtaPhi->GetXaxis()->SetTitle("phi"); h2_missedEtaPhi->GetYaxis()->SetTitle("eta");
  //h2_missedEtaPhi->Draw("COLZ");
  //c4->SaveAs(Form("%s/MissedEtaPhi_ele_%s.png", plotDir.c_str(), output_base.c_str()));


  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/output_ele_%s.root", output_base.c_str()),"recreate");

  denom->Write();
  denom_endcap->Write();

  num_15->Write();
  num_20->Write();
  num_30->Write();
  num_40->Write();
  num_50->Write();

  l1_7->Write();
  l1_15->Write();
  l1_21->Write();

  num_15_endcap->Write();
  num_20_endcap->Write();
  num_30_endcap->Write();
  num_40_endcap->Write();
  num_50_endcap->Write();

  l1_7_endcap->Write();
  l1_15_endcap->Write();
  l1_21_endcap->Write();

  r_15->Write();
  r_20->Write();
  r_30->Write();
  r_40->Write();
  r_50->Write();

  r_15_endcap->Write();
  r_20_endcap->Write();
  r_30_endcap->Write();
  r_40_endcap->Write();
  r_50_endcap->Write();

  wf->Close();

}
