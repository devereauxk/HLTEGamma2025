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

TH1D *h_hiBin = new TH1D("h_hiBin","h_hiBin",250,0,250);



// ==========================================================
// main: trigger turn on curves for HLT photons in MC
// ==========================================================

void triggerAnalysis_pho_mc(

    // 2025 QCDPhoton sample
    string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/PhotonQCD_pTMin20_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_QCDPhoton/251023_035328/0000/",
    string inputHLT = "", //"/afs/cern.ch/user/k/kdeverea/HLTClayton/CMSSW_15_1_0/src/workstation/HLT_emulation/scripts/PbPb/openHLTfiles/",
    string output_base = "MCQCDPhoton30_100",

    int nfiles = 200,
    float minHiBin = 60.0,
    float maxHiBin = 200.0
  ){

  std::cout << "running triggerAnalysis_pho_mc()" << std::endl;
  std::cout << "input forest directory = " << inputForest  << std::endl;
  std::cout << "input HLT directory    = " << inputHLT  << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  TChain *Tree_HLT_HIGEDPhoton20 = new TChain("hltobject/HLT_HIGEDPhoton20_v");
  TChain *HltTree_emulated = new TChain("hltanalysis/HltTree");
  
  std::cout<< "Adding input files...";

  // forest files
  if (nfiles == -1) {
    string this_forest_name = inputForest + "/HiForestMiniAOD_*.root";
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

  // hlt files
  /*
  if (nfiles == -1) {
    string this_hlt_name = inputHLT + "/*openHLT_Zee_*.root";
    Tree_HLT_HIGEDPhoton20     ->Add(this_hlt_name.c_str());
    HltTree_emulated        ->Add(this_hlt_name.c_str());
  }
  else {
    for(int i = 0; i < nfiles; i++) {
      string this_hlt_name = inputHLT + "/20251030_openHLT_Zee_" + to_string(i) + ".root";
      Tree_HLT_HIGEDPhoton20   ->Add(this_hlt_name.c_str());
      HltTree_emulated      ->Add(this_hlt_name.c_str());
    }
  }
  */
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
  Int_t           HLT_HIGEDPhoton20;
  Int_t           HLT_HIGEDPhoton30;
  Int_t           HLT_HIGEDPhoton40;
  Int_t           HLT_HIGEDPhoton50;
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

  HltTree->SetBranchStatus("HLT_HIGEDPhoton20_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton30_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton40_v16", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton50_v16", 1);

  HltTree->SetBranchAddress("Event", &forest_Event);
  HltTree->SetBranchAddress("LumiBlock", &forest_LumiBlock);
  HltTree->SetBranchAddress("Run", &forest_Run);

  HltTree->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7); 
  HltTree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15);
  HltTree->SetBranchAddress("L1_SingleEG21_BptxAND", &L1_SingleEG21);

  HltTree->SetBranchAddress("HLT_HIGEDPhoton20_v16", &HLT_HIGEDPhoton20);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton30_v16", &HLT_HIGEDPhoton30);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton40_v16", &HLT_HIGEDPhoton40);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton50_v16", &HLT_HIGEDPhoton50);

  // hlt emulated
  /*
  ULong64_t       hlt_Event;
  Int_t           hlt_LumiBlock, hlt_Run;

  HltTree_emulated->SetBranchStatus("*",0);     // disable all branches
  HltTree_emulated->SetBranchStatus("Event", 1);
  HltTree_emulated->SetBranchStatus("LumiBlock", 1);
  HltTree_emulated->SetBranchStatus("Run", 1);

  HltTree_emulated->SetBranchAddress("Event", &hlt_Event);
  HltTree_emulated->SetBranchAddress("LumiBlock", &hlt_LumiBlock);
  HltTree_emulated->SetBranchAddress("Run", &hlt_Run);

  // HLT_HIGEDPhoton20 from emulated HLT
  vector<float>           *Ele20Gsf_pt = nullptr;
  vector<float>           *Ele20Gsf_eta = nullptr;
  vector<float>           *Ele20Gsf_phi = nullptr;

  Tree_HLT_HIGEDPhoton20->SetBranchStatus("*",0);     // disable all branches
  Tree_HLT_HIGEDPhoton20->SetBranchStatus("pt", 1);
  Tree_HLT_HIGEDPhoton20->SetBranchStatus("eta", 1);
  Tree_HLT_HIGEDPhoton20->SetBranchStatus("phi", 1);

  Tree_HLT_HIGEDPhoton20->SetBranchAddress("pt", &Ele20Gsf_pt);
  Tree_HLT_HIGEDPhoton20->SetBranchAddress("eta", &Ele20Gsf_eta);
  Tree_HLT_HIGEDPhoton20->SetBranchAddress("phi", &Ele20Gsf_phi);
  */


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

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

  // mc branches
  EventTree->SetBranchStatus("pho_genMatchedIndex", 1);
  EventTree->SetBranchStatus("mcEt", 1);
  EventTree->SetBranchStatus("mcPID",1);
  EventTree->SetBranchStatus("mcEta",1);
  EventTree->SetBranchStatus("mcPhi",1);

  EventTree->SetBranchStatus("nPho",1);
  EventTree->SetBranchStatus("phoEt",1);
  EventTree->SetBranchStatus("phoEta",1);
  EventTree->SetBranchStatus("phoHoverE",1);
  EventTree->SetBranchStatus("phoPhi",1);

  EventTree->SetBranchStatus("phoSCEta",1);
  EventTree->SetBranchStatus("phoSigmaIEtaIEta_2012",1);
  EventTree->SetBranchStatus("mcCalIsoDR04",1);
  EventTree->SetBranchStatus("pfpIso3subUE",1);
  EventTree->SetBranchStatus("pfcIso3subUE",1);
  EventTree->SetBranchStatus("pfnIso3subUE",1);
  EventTree->SetBranchStatus("pho_swissCrx",1);
  EventTree->SetBranchStatus("pho_seedTime",1);
  EventTree->SetBranchStatus("pho_ecalClusterIsoR3",1);
  EventTree->SetBranchStatus("pho_hcalRechitIsoR3",1);
  EventTree->SetBranchStatus("pho_trackIsoR3PtCut20",1);

  EventTree->SetBranchAddress("pho_genMatchedIndex", &pho_genMatchedIndex);
  EventTree->SetBranchAddress("mcEt", &mcEt);
  EventTree->SetBranchAddress("mcPID", &mcPID);
  EventTree->SetBranchAddress("mcEta", &mcEta);
  EventTree->SetBranchAddress("mcPhi", &mcPhi);

  EventTree->SetBranchAddress("nPho", &nPho);
  EventTree->SetBranchAddress("phoEt", &phoEt);
  EventTree->SetBranchAddress("phoEta", &phoEta);
  EventTree->SetBranchAddress("phoHoverE", &phoHoverE);
  EventTree->SetBranchAddress("phoPhi", &phoPhi);

  EventTree->SetBranchAddress("phoSCEta", &phoSCEta);
  EventTree->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012);
  EventTree->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
  EventTree->SetBranchAddress("pfpIso3subUE", &pfpIso3subUE);
  EventTree->SetBranchAddress("pfcIso3subUE", &pfcIso3subUE);
  EventTree->SetBranchAddress("pfnIso3subUE", &pfnIso3subUE);
  EventTree->SetBranchAddress("pho_swissCrx", &pho_swissCrx);
  EventTree->SetBranchAddress("pho_seedTime", &pho_seedTime);
  EventTree->SetBranchAddress("pho_ecalClusterIsoR3", &pho_ecalClusterIsoR3);
  EventTree->SetBranchAddress("pho_hcalRechitIsoR3", &pho_hcalRechitIsoR3);
  EventTree->SetBranchAddress("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20);


  std::cout << "done" << std::endl;

  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop  =================

  // loop through reco objects
  int n_reco = HltTree-> GetEntries();
  int nphos = 0;
  int npass_leading = 0;
  int npass_genmatch = 0;
  int npass_trkrcut = 0;
  int npass_l1 = 0;
  int npass_trigger = 0;
  for (ULong64_t i_event = 0; i_event < HltTree->GetEntries(); ++i_event){

    if(i_event%(HltTree->GetEntries()/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    HiTree->GetEntry(i_event);
    HltTree->GetEntry(i_event);
    EventTree->GetEntry(i_event);
  	
    h_hiBin->Fill(hiBin);

    // event cuts
    if(fabs(vz)>15.0) continue;
    if(hiBin < minHiBin || hiBin >= maxHiBin) continue;

    nphos += nPho;

    // ========== skip duplicate reco events ===========
    // TODO

    // ============= find leading electron ==============
    // take highest reco pT as leading electron
    float maxPt = 0;
    int i_leading = -1;
    
    // loop over reco
    for(Int_t i_track = 0; i_track < nPho; i_track++){

      // find gen matched photons
      int i_genmatch = pho_genMatchedIndex->at(i_track);
      if (i_genmatch < 0) continue;
      if (abs(mcPID->at(i_genmatch)) != 22) continue;

      // check if matches to hltobject
      /*
      bool hlt_matched = false;
      for(Int_t i_hlt = 0; i_hlt < Ele30Gsf_pt->size(); i_hlt++) {
        // deltaR requirement
        float dEta = eleEta->at(i_track) - Ele30Gsf_eta->at(i_hlt);
        float dPhi = elePhi->at(i_track) - Ele30Gsf_phi->at(i_hlt);
        if (dPhi > TMath::Pi()) dPhi -= TMath::Pi();
        float dR = sqrt( dEta*dEta + dPhi*dPhi );

        if (dR < 0.5) {
          hlt_matched = true;
          break;
        }
      }
        */

      // compare with previous leading candidate
      if(phoEt->at(i_track) > maxPt) { // find the leading phoEt in the event
        maxPt = phoEt->at(i_track);
        i_leading = i_track;
      }

    }
    if (i_leading == -1) continue;
    npass_leading++;


    // ============= apply track rejection =============
    float sumIso = pfpIso3subUE->at(i_leading) + pfcIso3subUE->at(i_leading) + pfnIso3subUE->at(i_leading);

    bool isBarrel = false;
    bool isEndcap = false;

    // barrel cut
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

    npass_trkrcut++;


    // =============== fill histograms ================

    float weight = pthat_weight;

    if(isBarrel) {

      // fill denominator histograms
      denom->Fill(maxPt, weight);
      
      if(L1_SingleEG7) {
        l1_7->Fill(maxPt, weight);
        npass_l1++;
      }

      if(L1_SingleEG15) l1_15->Fill(maxPt, weight);
      if(L1_SingleEG21) l1_21->Fill(maxPt, weight);

      // fill numerator histograms
      if(HLT_HIGEDPhoton20) {
        num_20->Fill(maxPt, weight);
        npass_trigger++;
      }

      if(HLT_HIGEDPhoton30) num_30->Fill(maxPt, weight);

      if(HLT_HIGEDPhoton40) num_40->Fill(maxPt, weight);

      if(HLT_HIGEDPhoton50) num_50->Fill(maxPt, weight);
    
    }

    if(isEndcap) {

      // fill denominator histograms
      denom_endcap->Fill(maxPt, weight);
      
      if(L1_SingleEG7) {
        l1_7_endcap->Fill(maxPt, weight);
      }

      if(L1_SingleEG15) l1_15_endcap->Fill(maxPt, weight);
      if(L1_SingleEG21) l1_21_endcap->Fill(maxPt, weight);

      // fill numerator histograms
      if(HLT_HIGEDPhoton20) {
        num_20_endcap->Fill(maxPt, weight);
      }

      if(HLT_HIGEDPhoton30) num_30_endcap->Fill(maxPt, weight);

      if(HLT_HIGEDPhoton40) num_40_endcap->Fill(maxPt, weight);

      if(HLT_HIGEDPhoton50) num_50_endcap->Fill(maxPt, weight);
    
    }
    


  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIGEDPhoton30 results ... " << std::endl;
  std::cout << "nphos = " << nphos << std::endl;
  std::cout << "npass_leading = " << npass_leading  << std::endl;
  std::cout << "npass_trkrcut = " << npass_trkrcut  << std::endl;
  std::cout << "npass_trigger = " << npass_trigger  << std::endl;

  /*
  r_20->Divide(num_20,denom,1,1,"B");
  r_30->Divide(num_30,denom,1,1,"B");
  r_40->Divide(num_40,denom,1,1,"B");
  r_50->Divide(num_50,denom,1,1,"B");
  */

  r_20->Divide(num_20,l1_7,1,1,"B");
  r_30->Divide(num_30,l1_15,1,1,"B");
  r_40->Divide(num_40,l1_21,1,1,"B");
  r_50->Divide(num_50,l1_21,1,1,"B"); 

  rl1_7->Divide(l1_7,denom,1,1,"B");
  rl1_15->Divide(l1_15,denom,1,1,"B");
  rl1_21->Divide(l1_21,denom,1,1,"B");

  r_20_endcap->Divide(num_20_endcap,l1_7_endcap,1,1,"B");
  r_30_endcap->Divide(num_30_endcap,l1_15_endcap,1,1,"B");
  r_40_endcap->Divide(num_40_endcap,l1_21_endcap,1,1,"B");
  r_50_endcap->Divide(num_50_endcap,l1_21_endcap,1,1,"B");

  r_20->SetLineColor(kRed-4);
  r_30->SetLineColor(kBlue-4);
  r_40->SetLineColor(kGreen+2);
  r_50->SetLineColor(kMagenta-9);
  r_20_endcap->SetLineColor(kRed-4);
  r_30_endcap->SetLineColor(kBlue-4);
  r_40_endcap->SetLineColor(kGreen+2);
  r_50_endcap->SetLineColor(kMagenta-9);

  r_20->SetMarkerColor(kRed-4);
  r_30->SetMarkerColor(kBlue-4);
  r_40->SetMarkerColor(kGreen+2);
  r_50->SetMarkerColor(kMagenta-9);
  r_20_endcap->SetMarkerColor(kRed-4);
  r_30_endcap->SetMarkerColor(kBlue-4);
  r_40_endcap->SetMarkerColor(kGreen+2);
  r_50_endcap->SetMarkerColor(kMagenta-9);

  double line_width = 1.8;
  r_20->SetLineWidth(line_width);
  r_30->SetLineWidth(line_width);
  r_40->SetLineWidth(line_width);
  r_50->SetLineWidth(line_width);
  r_20_endcap->SetLineWidth(line_width);
  r_30_endcap->SetLineWidth(line_width);
  r_40_endcap->SetLineWidth(line_width);
  r_50_endcap->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_20->SetMarkerSize(marker_size);
  r_30->SetMarkerSize(marker_size);
  r_40->SetMarkerSize(marker_size);
  r_50->SetMarkerSize(marker_size);
  r_20_endcap->SetMarkerSize(marker_size);
  r_30_endcap->SetMarkerSize(marker_size);
  r_40_endcap->SetMarkerSize(marker_size);
  r_50_endcap->SetMarkerSize(marker_size);
    
  r_20->SetMarkerStyle(20);
  r_30->SetMarkerStyle(21);
  r_40->SetMarkerStyle(22);
  r_50->SetMarkerStyle(23);
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
  r_20->GetXaxis()->SetTitle("photon #font[52]{E}_{T} [GeV]");
  r_20->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_20,"HLT_HIGEDPhoton20_v16");
  leg->AddEntry(r_30,"HLT_HIGEDPhoton30_v16");
  leg->AddEntry(r_40,"HLT_HIGEDPhoton40_v16");
  leg->AddEntry(r_50,"HLT_HIGEDPhoton50_v16");
  //leg->SetBorderSize(0);
  r_20->Draw();
  leg->Draw();
  r_30->Draw("same");
  r_40->Draw("same");
  r_50->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet QCDPhoton");
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("plots/TriggerEfficiency_pho_%s_barrel.png", output_base.c_str()));


  // =============== Fig HLT (no L1 emulation) ENDCAP ==================
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
  r_20_endcap->GetXaxis()->SetTitle("photon #font[52]{E}_{T} [GeV]");
  r_20_endcap->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg_endcap = new TLegend(0.55,0.3,0.88,0.5);
  leg_endcap->AddEntry(r_20_endcap,"HLT_HIGEDPhoton20_v16");
  leg_endcap->AddEntry(r_30_endcap,"HLT_HIGEDPhoton30_v16");
  leg_endcap->AddEntry(r_40_endcap,"HLT_HIGEDPhoton40_v16");
  leg_endcap->AddEntry(r_50_endcap,"HLT_HIGEDPhoton50_v16");
  //leg_endcap->SetBorderSize(0);
  r_20_endcap->Draw();
  leg_endcap->Draw();
  r_30_endcap->Draw("same");
  r_40_endcap->Draw("same");
  r_50_endcap->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet QCDPhoton");
  la_endcap->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la_endcap->DrawLatexNDC(0.6,0.65,"Endcap: 1.556 < |#eta| < 2.1");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("plots/TriggerEfficiency_pho_%s_endcap.png", output_base.c_str()));


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

  la2->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet QCDPhoton");
  la2->DrawLatexNDC(0.72,0.92,"L1 emulation");
  la2->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la2->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c2->SaveAs(Form("plots/TriggerEfficiency_pho_%s_L1.png", output_base.c_str()));


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
  c3->SaveAs(Form("plots/hiBin_pho_%s.png", output_base.c_str()));


  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/output_pho_%s.root", output_base.c_str()),"recreate");

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

  r_20->Write();
  r_30->Write();
  r_40->Write();
  r_50->Write();

  r_20_endcap->Write();
  r_30_endcap->Write();
  r_40_endcap->Write();
  r_50_endcap->Write();

  wf->Close();

}


