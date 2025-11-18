using namespace std;

Float_t PtMin = 2.;
Float_t PtMax = 110;
Int_t NPtBins = 50;

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

TH1D *ele = new TH1D("ele","ele",NPtBins,PtMin,PtMax);

TH1D *num_10 = new TH1D("num_10","num_10",NPtBins,PtMin,PtMax);
TH1D *num_10_endcap = new TH1D("num_10_endcap","num_10_endcap",NPtBins,PtMin,PtMax);
TH1D *r_10 = new TH1D("r_10","r_10",NPtBins,PtMin,PtMax);
TH1D *r_10_endcap = new TH1D("r_10_endcap","r_10_endcap",NPtBins,PtMin,PtMax);

TH1D *num_20 = new TH1D("num_20","num_20",NPtBins,PtMin,PtMax);
TH1D *num_20_endcap = new TH1D("num_20_endcap","num_20_endcap",NPtBins,PtMin,PtMax);
TH1D *r_pho10_ele = new TH1D("r_pho10_ele","r_pho10_ele",NPtBins,PtMin,PtMax);
TH1D *r_pho10_ele_endcap = new TH1D("r_pho10_ele_endcap","r_pho10_ele_endcap",NPtBins,PtMin,PtMax);

TH1D *num_30 = new TH1D("num_30","num_30",NPtBins,PtMin,PtMax);
TH1D *num_30_endcap = new TH1D("num_30_endcap","num_30_endcap",NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele = new TH1D("r_pho20_ele","r_pho20_ele",NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele_endcap = new TH1D("r_pho20_ele_endcap","r_pho20_ele_endcap",NPtBins,PtMin,PtMax);

TH1D *num_40 = new TH1D("num_40","num_40",NPtBins,PtMin,PtMax);
TH1D *num_40_endcap = new TH1D("num_40_endcap","num_40_endcap",NPtBins,PtMin,PtMax);
TH1D *r_pho30_ele = new TH1D("r_pho30_ele","r_pho30_ele",NPtBins,PtMin,PtMax);
TH1D *r_pho30_ele_endcap = new TH1D("r_pho30_ele_endcap","r_pho30_ele_endcap",NPtBins,PtMin,PtMax);

TH1D *num_50 = new TH1D("num_50","num_50",NPtBins,PtMin,PtMax);
TH1D *num_50_endcap = new TH1D("num_50_endcap","num_50_endcap",NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele20 = new TH1D("r_pho20_ele20","r_pho20_ele20",NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele20_endcap = new TH1D("r_pho20_ele20_endcap","r_pho20_ele20_endcap",NPtBins,PtMin,PtMax);

TH1D *h_hiBin = new TH1D("h_hiBin","h_hiBin",250,0,250);



// ==========================================================
// main: trigger turn on curves for HLT photons in PbPb MC
// ==========================================================

void triggerAnalysis_pho_over_ele_PbPb(

    // PbPb data
    string inputForest = "/eos/cms/store/group/phys_heavyions/nbarnett/Forests/run399499/HIPhysicsRawPrime0/CRAB_UserFiles/crab_forest_PbPb_RP0_run399466_11_15_2025_v1/251115_200501/0000",
    string inputText = "run399499_forests.txt",
    string output_base = "DataPbPb0_100_run399499",
    string plot_label = "Run 399499, HIPhysicsRawPrime0-14",

    int nfiles = -1,
    float minHiBin = 0.0,
    float maxHiBin = 200.0
  ){

  std::cout << "running triggerAnalysis_pho_PbPb()" << std::endl;
  std::cout << "input forest directory = " << inputForest  << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");
  
  std::cout<< "Adding input files...";

  // forest files
  if (inputText != "") {
    ifstream infile(inputText);
    string line;
    int fileCount = 0;
    while (getline(infile, line) && (fileCount < nfiles || nfiles == -1)) {
      HltTree   ->Add(line.c_str());
      EventTree ->Add(line.c_str());
      HiTree    ->Add(line.c_str());
      fileCount++;
    }
  }
  else if (nfiles == -1) {
    string this_forest_name = inputForest + "/HiForest*.root";
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
  HiTree->SetBranchStatus("vz",1);

  HiTree->SetBranchAddress("evt", &evt);
  HiTree->SetBranchAddress("lumi", &lumi);
  HiTree->SetBranchAddress("run", &run);
  HiTree->SetBranchAddress("hiBin", &hiBin);
  HiTree->SetBranchAddress("hiHF", &hiHF);
  HiTree->SetBranchAddress("vz", &vz);


  // ================= HLT Trees =================
  // forest
  Int_t           L1_SingleEG7;
  Int_t           L1_SingleEG15;
  Int_t           L1_SingleEG21;
  Int_t           HLT_HIMinimumBias;
  Int_t           HLT_HIGEDPhoton10;
  Int_t           HLT_HIGEDPhoton20;
  Int_t           HLT_HIGEDPhoton30;
  Int_t           HLT_HIGEDPhoton40;
  Int_t           HLT_HIGEDPhoton50;
  Int_t           HLT_HIEle10Gsf;
  Int_t           HLT_HIEle15Gsf;
  Int_t           HLT_HIEle20Gsf;
  Int_t           HLT_HIEle30Gsf;
  Int_t           HLT_HIEle40Gsf;
  Int_t           HLT_HIEle50Gsf;
  ULong64_t       forest_Event;
  Int_t           forest_LumiBlock, forest_Run;

  //HltTree->Print();
  //HltTree->Show(0);

  HltTree->SetBranchStatus("*",0);     // disable all branches
  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  HltTree->SetBranchStatus("L1_SingleEG7_BptxAND", 1);
  HltTree->SetBranchStatus("L1_SingleEG15_BptxAND", 1);
  HltTree->SetBranchStatus("L1_SingleEG21_BptxAND", 1);

  HltTree->SetBranchStatus("HLT_HIMinimumBiasHF1ANDZDC1nOR_v6", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton10_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton20_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton30_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton40_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton50_v16", 1);

  HltTree->SetBranchStatus("HLT_HIEle10Gsf_v17", 1);
  HltTree->SetBranchStatus("HLT_HIEle15Gsf_v17", 1);
  HltTree->SetBranchStatus("HLT_HIEle20Gsf_v17", 1);
  HltTree->SetBranchStatus("HLT_HIEle30Gsf_v17", 1);
  HltTree->SetBranchStatus("HLT_HIEle40Gsf_v17", 1);
  HltTree->SetBranchStatus("HLT_HIEle50Gsf_v17", 1);

  HltTree->SetBranchAddress("Event", &forest_Event);
  HltTree->SetBranchAddress("LumiBlock", &forest_LumiBlock);
  HltTree->SetBranchAddress("Run", &forest_Run);

  HltTree->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7); 
  HltTree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15);
  HltTree->SetBranchAddress("L1_SingleEG21_BptxAND", &L1_SingleEG21);

  HltTree->SetBranchAddress("HLT_HIMinimumBiasHF1ANDZDC1nOR_v6", &HLT_HIMinimumBias);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton10_v16", &HLT_HIGEDPhoton10);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton20_v16", &HLT_HIGEDPhoton20);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton30_v16", &HLT_HIGEDPhoton30);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton40_v16", &HLT_HIGEDPhoton40);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton50_v16", &HLT_HIGEDPhoton50);

  HltTree->SetBranchAddress("HLT_HIEle10Gsf_v17", &HLT_HIEle10Gsf);
  HltTree->SetBranchAddress("HLT_HIEle15Gsf_v17", &HLT_HIEle15Gsf);
  HltTree->SetBranchAddress("HLT_HIEle20Gsf_v17", &HLT_HIEle20Gsf);
  HltTree->SetBranchAddress("HLT_HIEle30Gsf_v17", &HLT_HIEle30Gsf);
  HltTree->SetBranchAddress("HLT_HIEle40Gsf_v17", &HLT_HIEle40Gsf);
  HltTree->SetBranchAddress("HLT_HIEle50Gsf_v17", &HLT_HIEle50Gsf);


  // ================= Event Tree =================
  vector<float>   *phoEt = nullptr;
  vector<float>   *phoEta = nullptr;
  vector<float>   *phoPhi = nullptr;
  Int_t           nPho;
  vector<int>     *pho_genMatchedIndex = nullptr;
  vector<float>   *mcEt = nullptr;
  vector<int>     *mcPID = nullptr;
  vector<float>   *mcEta = nullptr;
  vector<float>   *mcPhi = nullptr;
  vector<float>   *phoHoverE = nullptr;
  vector<float>   *phoSCEta = nullptr;
  vector<float>   *phoSigmaIEtaIEta_2012 = nullptr;
  vector<float>   *mcCalIsoDR04 = nullptr;
  vector<float>   *pfpIso3subUE = nullptr;
  vector<float>   *pfcIso3subUE = nullptr;
  vector<float>   *pfnIso3subUE = nullptr;
  vector<float>   *pho_swissCrx = nullptr;
  vector<float>   *pho_seedTime = nullptr;
  vector<float>   *pho_ecalClusterIsoR3 = nullptr;
  vector<float>   *pho_hcalRechitIsoR3 = nullptr;
  vector<float>   *pho_trackIsoR3PtCut20 = nullptr;

  vector<float>   *elePt = nullptr;
  vector<float>   *eleEta = nullptr;
  vector<float>   *elePhi = nullptr;
  int             nEle;
  vector<int>     *ele_genMatchedIndex = nullptr;
  vector<float>   *eleHoverE = nullptr;
  vector<float>   *eleSigmaIEtaIEta_2012 = nullptr;
  vector<float>   *eledEtaSeedAtVtx = nullptr;
  vector<float>   *eledPhiAtVtx = nullptr;
  vector<float>   *eleEoverPInv = nullptr;
  vector<int>     *eleMissHits = nullptr;
  vector<float>   *eleIP3D = nullptr;
  vector<float>   *eleSCEta = nullptr;
  vector<float>   *eleSCRawEn = nullptr;
  

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

  // mc branches

  EventTree->SetBranchStatus("nPho",1);
  EventTree->SetBranchStatus("phoEt",1);
  EventTree->SetBranchStatus("phoEta",1);
  EventTree->SetBranchStatus("phoHoverE",1);
  EventTree->SetBranchStatus("phoPhi",1);

  EventTree->SetBranchStatus("phoSCEta",1);
  EventTree->SetBranchStatus("phoSigmaIEtaIEta_2012",1);
  //EventTree->SetBranchStatus("mcCalIsoDR04",1);
  EventTree->SetBranchStatus("pfpIso3subUE",1);
  EventTree->SetBranchStatus("pfcIso3subUE",1);
  EventTree->SetBranchStatus("pfnIso3subUE",1);
  EventTree->SetBranchStatus("pho_swissCrx",1);
  EventTree->SetBranchStatus("pho_seedTime",1);
  EventTree->SetBranchStatus("pho_ecalClusterIsoR3",1);
  EventTree->SetBranchStatus("pho_hcalRechitIsoR3",1);
  EventTree->SetBranchStatus("pho_trackIsoR3PtCut20",1);

  EventTree->SetBranchAddress("nPho", &nPho);
  EventTree->SetBranchAddress("phoEt", &phoEt);
  EventTree->SetBranchAddress("phoEta", &phoEta);
  EventTree->SetBranchAddress("phoHoverE", &phoHoverE);
  EventTree->SetBranchAddress("phoPhi", &phoPhi);

  EventTree->SetBranchAddress("phoSCEta", &phoSCEta);
  EventTree->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012);
  //EventTree->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
  EventTree->SetBranchAddress("pfpIso3subUE", &pfpIso3subUE);
  EventTree->SetBranchAddress("pfcIso3subUE", &pfcIso3subUE);
  EventTree->SetBranchAddress("pfnIso3subUE", &pfnIso3subUE);
  EventTree->SetBranchAddress("pho_swissCrx", &pho_swissCrx);
  EventTree->SetBranchAddress("pho_seedTime", &pho_seedTime);
  EventTree->SetBranchAddress("pho_ecalClusterIsoR3", &pho_ecalClusterIsoR3);
  EventTree->SetBranchAddress("pho_hcalRechitIsoR3", &pho_hcalRechitIsoR3);
  EventTree->SetBranchAddress("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20);

  EventTree->SetBranchStatus("nEle",1);
  EventTree->SetBranchStatus("elePt",1);
  EventTree->SetBranchStatus("eleEta",1);
  EventTree->SetBranchStatus("eleHoverE",1);
  EventTree->SetBranchStatus("elePhi",1);

  EventTree->SetBranchStatus("eleSigmaIEtaIEta_2012",1);
  EventTree->SetBranchStatus("eledEtaSeedAtVtx",1);
  EventTree->SetBranchStatus("eledPhiAtVtx",1);
  EventTree->SetBranchStatus("eleEoverPInv",1);
  EventTree->SetBranchStatus("eleMissHits",1);
  EventTree->SetBranchStatus("eleIP3D",1);
  EventTree->SetBranchStatus("eleSCEta",1);
  EventTree->SetBranchStatus("eleSCRawEn",1);

  EventTree->SetBranchAddress("nEle", &nEle);
  EventTree->SetBranchAddress("elePt", &elePt);
  EventTree->SetBranchAddress("eleEta", &eleEta);
  EventTree->SetBranchAddress("eleHoverE", &eleHoverE);
  EventTree->SetBranchAddress("elePhi", &elePhi);

  EventTree->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
  EventTree->SetBranchAddress("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
  EventTree->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
  EventTree->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
  EventTree->SetBranchAddress("eleMissHits", &eleMissHits);
  EventTree->SetBranchAddress("eleIP3D", &eleIP3D);
  EventTree->SetBranchAddress("eleSCEta", &eleSCEta);
  EventTree->SetBranchAddress("eleSCRawEn", &eleSCRawEn);


  std::cout << "done" << std::endl;
  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop  =================

  // loop through reco objects
  int n_reco = HltTree-> GetEntries();
  int neles = 0;
  int npass_leading = 0;
  int npass_trkrcut = 0;
  int npass_denom = 0;
  int npass_numer = 0;
  for (ULong64_t i_event = 0; i_event < HltTree->GetEntries(); ++i_event){

    if(i_event%(HltTree->GetEntries()/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    HiTree->GetEntry(i_event);
    HltTree->GetEntry(i_event);
    EventTree->GetEntry(i_event);
  	
    h_hiBin->Fill(hiBin);

    // event cuts
    if(fabs(vz)>15.0) continue;
    if(hiBin < minHiBin || hiBin >= maxHiBin) continue;

    neles += nEle;

    // ========== skip duplicate reco events ===========
    // NO

    // ===== find HLT emulation match of reco event =====
    // NO


    // ============= find leading photon ==============
    // take highest reco pT as leading electron
    /*
    float maxPt = 0;
    int i_leading = -1;
    
    // loop over reco
    for(Int_t i_track = 0; i_track < nPho; i_track++){

      // compare with previous leading candidate
      if(phoEt->at(i_track) > maxPt) { // find the leading phoEt in the event
        maxPt = phoEt->at(i_track);
        i_leading = i_track;
      }

    }
    if (i_leading == -1) continue;
    npass_leading++;
    */


    // =============== find leading electron ==============
    float maxPt = 0;
    int i_leading = -1;

    // loop over reco
    for(Int_t i_ele = 0; i_ele < nEle; i_ele++){

      // compare with previous leading candidate
      if(elePt->at(i_ele) > maxPt) { // find the leading elePt in the event
        maxPt = elePt->at(i_ele);
        i_leading = i_ele;
      }

    }
    if (i_leading == -1) continue;
    npass_leading++;


    // ============= apply track rejection =============
    //float sumIso = pfpIso3subUE->at(i_leading) + pfcIso3subUE->at(i_leading) + pfnIso3subUE->at(i_leading);

    bool isBarrel_ele = false;

    // barrel cut
    /*
    if (abs(phoSCEta->at(i_leading)) < 1.442) {
      if (sumIso > 5)                                      continue;
      if (phoHoverE->at(i_leading) > 0.2)                  continue;
      if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.012)    continue;
      if (pho_swissCrx->at(i_leading) > 0.9)               continue;
      if (fabs(pho_seedTime->at(i_leading)) > 3.0)         continue;
      isBarrel = true;
    }

    // endcap cut
    else if (abs(phoSCEta->at(i_leading)) > 1.556 && abs(phoSCEta->at(i_leading)) < 2.1) {
      if (sumIso > 10)                                     continue;
      if (phoHoverE->at(i_leading) > 0.3)                  continue;
      if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.03)     continue;
      if (pho_swissCrx->at(i_leading) > 0.9)               continue;
      if (fabs(pho_seedTime->at(i_leading)) > 3.0)         continue;
      isEndcap = true;
    }
    else continue;
    */

    if (abs(eleEta->at(i_leading)) < 1.442) {
      if (eleMissHits->at(i_leading) > 3)                 continue;
      if (eleIP3D->at(i_leading) > 0.03)                  continue;
      if (eleSigmaIEtaIEta_2012->at(i_leading) > 0.012)   continue;
      if (eledEtaSeedAtVtx->at(i_leading) > 0.0037)       continue;
      if (eledPhiAtVtx->at(i_leading) > 0.1280)           continue;
      if (eleEoverPInv->at(i_leading) > 0.1065)           continue;
      if (eleHoverE->at(i_leading) > 0.13)                continue;
      isBarrel_ele = true;
    }
    else continue;

    npass_trkrcut++;


    // =============== match to gen level ===============
    // NO


    // =============== fill histograms ================

    float weight = 1;

    if(isBarrel_ele) {

      // fill denominator histograms
      denom->Fill(maxPt, weight);
      
      if(L1_SingleEG7) l1_7->Fill(maxPt, weight);
      if(L1_SingleEG15) l1_15->Fill(maxPt, weight);
      if(L1_SingleEG21) l1_21->Fill(maxPt, weight);


      // fill numerator (photon) histograms
      if(HLT_HIMinimumBias && HLT_HIGEDPhoton10) {
        num_10->Fill(maxPt, weight);
        npass_numer++;
      }

      if(HLT_HIMinimumBias && HLT_HIGEDPhoton20) num_20->Fill(maxPt, weight);

      if(HLT_HIMinimumBias && HLT_HIGEDPhoton30) num_30->Fill(maxPt, weight);

      if(HLT_HIMinimumBias && HLT_HIGEDPhoton40) num_40->Fill(maxPt, weight);

      if(HLT_HIMinimumBias && HLT_HIGEDPhoton50) num_50->Fill(maxPt, weight);


      // fill denominator (electron) histograms
      if(HLT_HIMinimumBias) {
        ele->Fill(maxPt, weight);
        npass_denom++;
      }
    
    }


  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIGEDPhoton10 results ... " << std::endl;
  std::cout << "neles = " << neles << std::endl;
  std::cout << "npass_leading = " << npass_leading  << std::endl;
  std::cout << "npass_trkrcut = " << npass_trkrcut  << std::endl;
  std::cout << "npass_denom = " << npass_denom  << std::endl;
  std::cout << "npass_numer = " << npass_numer  << std::endl;

  r_pho10_ele->Divide(num_10,ele,1,1,"B");
  r_pho20_ele->Divide(num_20,ele,1,1,"B");
  r_pho30_ele->Divide(num_30,ele,1,1,"B");

  r_pho10_ele->SetLineColor(kRed-4);
  r_pho20_ele->SetLineColor(kBlue-4);
  r_pho30_ele->SetLineColor(kGreen+2);
  r_pho20_ele20->SetLineColor(kMagenta-9);
  r_pho10_ele_endcap->SetLineColor(kRed-4);
  r_pho20_ele_endcap->SetLineColor(kBlue-4);
  r_pho30_ele_endcap->SetLineColor(kGreen+2);
  r_pho20_ele20_endcap->SetLineColor(kMagenta-9);

  r_pho10_ele->SetMarkerColor(kRed-4);
  r_pho20_ele->SetMarkerColor(kBlue-4);
  r_pho30_ele->SetMarkerColor(kGreen+2);
  r_pho20_ele20->SetMarkerColor(kMagenta-9);
  r_pho10_ele_endcap->SetMarkerColor(kRed-4);
  r_pho20_ele_endcap->SetMarkerColor(kBlue-4);
  r_pho30_ele_endcap->SetMarkerColor(kGreen+2);
  r_pho20_ele20_endcap->SetMarkerColor(kMagenta-9);

  double line_width = 1.8;
  r_pho10_ele->SetLineWidth(line_width);
  r_pho20_ele->SetLineWidth(line_width);
  r_pho30_ele->SetLineWidth(line_width);
  r_pho20_ele20->SetLineWidth(line_width);
  r_pho10_ele_endcap->SetLineWidth(line_width);
  r_pho20_ele_endcap->SetLineWidth(line_width);
  r_pho30_ele_endcap->SetLineWidth(line_width);
  r_pho20_ele20_endcap->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_pho10_ele->SetMarkerSize(marker_size);
  r_pho20_ele->SetMarkerSize(marker_size);
  r_pho30_ele->SetMarkerSize(marker_size);
  r_pho20_ele20->SetMarkerSize(marker_size);
  r_pho10_ele_endcap->SetMarkerSize(marker_size);
  r_pho20_ele_endcap->SetMarkerSize(marker_size);
  r_pho30_ele_endcap->SetMarkerSize(marker_size);
  r_pho20_ele20_endcap->SetMarkerSize(marker_size);
    
  r_pho10_ele->SetMarkerStyle(20);
  r_pho20_ele->SetMarkerStyle(21);
  r_pho30_ele->SetMarkerStyle(22);
  r_pho20_ele20->SetMarkerStyle(23);
  r_pho10_ele_endcap->SetMarkerStyle(20);
  r_pho20_ele_endcap->SetMarkerStyle(21);
  r_pho30_ele_endcap->SetMarkerStyle(22);
  r_pho20_ele20_endcap->SetMarkerStyle(23);

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
  r_pho10_ele->SetTitle("");
  r_pho10_ele->SetStats(0);
  r_pho10_ele->GetXaxis()->SetTitleSize(0.05);
  r_pho10_ele->GetYaxis()->SetTitleSize(0.05);
  r_pho10_ele->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_pho10_ele->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_pho10_ele,"HLT_HIGEDPhoton10");
  leg->AddEntry(r_pho20_ele,"HLT_HIGEDPhoton20");
  leg->AddEntry(r_pho30_ele,"HLT_HIGEDPhoton30");
  //leg->SetBorderSize(0);
  r_pho10_ele->Draw();
  leg->Draw();
  r_pho20_ele->Draw("same");
  r_pho30_ele->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,plot_label.c_str());
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("plots/TriggerEfficiency_phoOverEle_%s_barrel.png", output_base.c_str()));


  // =============== Fig HLT (no L1 emulation) ENDCAP ==================
  /*
  TCanvas *c1_endcap = new TCanvas("c1_endcap","c1_endcap",700,600);
  c1_endcap->cd();
  TPad *p1_endcap = new TPad("p1_endcap","p1_endcap",0,0,1,1);
  p1_endcap->SetLeftMargin(0.13);
  p1_endcap->SetBottomMargin(0.14);
  p1_endcap->Draw();
  p1_endcap->cd();
  r_pho10_ele_endcap->SetTitle("");
  r_pho10_ele_endcap->SetStats(0);
  r_pho10_ele_endcap->GetXaxis()->SetTitleSize(0.05);
  r_pho10_ele_endcap->GetYaxis()->SetTitleSize(0.05);
  r_pho10_ele_endcap->GetXaxis()->SetTitle("photon #font[52]{E}_{T} [GeV]");
  r_pho10_ele_endcap->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg_endcap = new TLegend(0.55,0.3,0.88,0.5);
  leg_endcap->AddEntry(r_pho10_ele_endcap,"HLT_HIGEDPhoton20_v16");
  leg_endcap->AddEntry(r_pho20_ele_endcap,"HLT_HIGEDPhoton30_v16");
  leg_endcap->AddEntry(r_pho30_ele_endcap,"HLT_HIGEDPhoton40_v16");
  leg_endcap->AddEntry(r_pho20_ele20_endcap,"HLT_HIGEDPhoton50_v16");
  //leg_endcap->SetBorderSize(0);
  r_pho10_ele_endcap->Draw();
  leg_endcap->Draw();
  r_pho20_ele_endcap->Draw("same");
  r_pho30_ele_endcap->Draw("same");
  r_pho20_ele20_endcap->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,plot_label.c_str());
  la_endcap->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la_endcap->DrawLatexNDC(0.6,0.65,"Endcap: 1.556 < |#eta| < 2.1");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("plots/TriggerEfficiency_phoOverEle_%s_endcap.png", output_base.c_str()));


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
  rl1_7->GetXaxis()->SetTitle("photon #font[52]{E}_{T} [GeV]");
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

  la2->DrawLatexNDC(0.22,0.92,plot_label.c_str());
  la2->DrawLatexNDC(0.72,0.92,"L1 emulation");
  la2->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la2->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c2->SaveAs(Form("plots/TriggerEfficiency_phoOverEle_%s_L1.png", output_base.c_str()));


  // ====== hiBin historogram simple ======
  TCanvas *c3 = new TCanvas("c3","c3",700,600);
  c3->cd();
  TPad *p3 = new TPad("p3","p3",0,0,1,1);
  p3->SetLeftMargin(0.13);
  p3->SetBottomMargin(0.14);
  p3->Draw();
  p3->cd();
  h_hiBin->SetTitle("");
  h_hiBin->SetStats(0);
  h_hiBin->GetXaxis()->SetTitleSize(0.05);
  h_hiBin->GetYaxis()->SetTitleSize(0.05);
  h_hiBin->GetXaxis()->SetTitle("hiBin");
  h_hiBin->GetYaxis()->SetTitle("Counts");
  h_hiBin->Draw();
  c3->SaveAs(Form("plots/hiBin_phoOverEle_%s.png", output_base.c_str()));
  */

  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/output_phoOverEle_%s.root", output_base.c_str()),"recreate");

  denom->Write();
  denom_endcap->Write();

  num_20->Write();
  num_30->Write();
  num_40->Write();
  num_50->Write();

  l1_7->Write();
  l1_15->Write();
  l1_21->Write();

  num_20_endcap->Write();
  num_30_endcap->Write();
  num_40_endcap->Write();
  num_50_endcap->Write();

  l1_7_endcap->Write();
  l1_15_endcap->Write();
  l1_21_endcap->Write();

  r_pho10_ele->Write();
  r_pho20_ele->Write();
  r_pho30_ele->Write();
  r_pho20_ele20->Write();

  r_pho10_ele_endcap->Write();
  r_pho20_ele_endcap->Write();
  r_pho30_ele_endcap->Write();
  r_pho20_ele20_endcap->Write();

  wf->Close();

}


