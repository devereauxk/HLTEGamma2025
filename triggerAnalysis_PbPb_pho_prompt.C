#include "include/EventMatcher.h"

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
TH1D *denomEta_40 = new TH1D("denomEta_40","denomEta_40",NEtaBins,etaMin,etaMax);
TH1D *denomEta_60 = new TH1D("denomEta_60","denomEta_60",NEtaBins,etaMin,etaMax);
TH1D *denomEta_70 = new TH1D("denomEta_70","denomEta_70",NEtaBins,etaMin,etaMax);
TH1D *denomEta_80 = new TH1D("denomEta_80","denomEta_80",NEtaBins,etaMin,etaMax);
TH1D *denomEta_100 = new TH1D("denomEta_100","denomEta_100",NEtaBins,etaMin,etaMax);
TH1D *denomEta_120 = new TH1D("denomEta_120","denomEta_120",NEtaBins,etaMin,etaMax);
TH1D *denomPhi_40 = new TH1D("denomPhi_40","denomPhi_40",NPhiBins,phiMin,phiMax);
TH1D *denomPhi_60 = new TH1D("denomPhi_60","denomPhi_60",NPhiBins,phiMin,phiMax);
TH1D *denomPhi_70 = new TH1D("denomPhi_70","denomPhi_70",NPhiBins,phiMin,phiMax);
TH1D *denomPhi_80 = new TH1D("denomPhi_80","denomPhi_80",NPhiBins,phiMin,phiMax);
TH1D *denomPhi_100 = new TH1D("denomPhi_100","denomPhi_100",NPhiBins,phiMin,phiMax);
TH1D *denomPhi_120 = new TH1D("denomPhi_120","denomPhi_120",NPhiBins,phiMin,phiMax);

TH1D *l1_7 = new TH1D("l1_7","l1_7",NPtBins,PtMin,PtMax);
TH1D *l1_15 = new TH1D("l1_15","l1_15",NPtBins,PtMin,PtMax);
TH1D *l1_21 = new TH1D("l1_21","l1_21",NPtBins,PtMin,PtMax);

TH1D *rl1_7 = new TH1D("rl1_7","rl1_7",NPtBins,PtMin,PtMax);
TH1D *rl1_15 = new TH1D("rl1_15","rl1_15",NPtBins,PtMin,PtMax);
TH1D *rl1_21 = new TH1D("rl1_21","rl1_21",NPtBins,PtMin,PtMax);

TH1D *num_20 = new TH1D("num_20","num_20",NPtBins,PtMin,PtMax);
TH1D *numEta_40 = new TH1D("numEta_40","numEta_40",NEtaBins,etaMin,etaMax);
TH1D *numPhi_40 = new TH1D("numPhi_40","numPhi_40",NPhiBins,phiMin,phiMax);
TH1D *r_20 = new TH1D("r_20","r_20",NPtBins,PtMin,PtMax);

TH1D *num_30 = new TH1D("num_30","num_30",NPtBins,PtMin,PtMax);
TH1D *numEta_60 = new TH1D("numEta_60","numEta_60",NEtaBins,etaMin,etaMax);
TH1D *numPhi_60 = new TH1D("numPhi_60","numPhi_60",NPhiBins,phiMin,phiMax);
TH1D *r_30 = new TH1D("r_30","r_30",NPtBins,PtMin,PtMax);

TH1D *num_40 = new TH1D("num_40","num_40",NPtBins,PtMin,PtMax);
TH1D *numEta_70 = new TH1D("numEta_70","numEta_70",NEtaBins,etaMin,etaMax);
TH1D *numPhi_70 = new TH1D("numPhi_70","numPhi_70",NPhiBins,phiMin,phiMax);
TH1D *r_40 = new TH1D("r_40","r_40",NPtBins,PtMin,PtMax);

TH1D *num_50 = new TH1D("num_50","num_50",NPtBins,PtMin,PtMax);
TH1D *numEta_100 = new TH1D("numEta_100","numEta_100",NEtaBins,etaMin,etaMax);
TH1D *numPhi_100 = new TH1D("numPhi_100","numPhi_100",NPhiBins,phiMin,phiMax);
TH1D *r_50 = new TH1D("r_50","r_50",NPtBins,PtMin,PtMax);

TH1D *num_120 = new TH1D("num_120","num_120",NPtBins,PtMin,PtMax);
TH1D *numEta_120 = new TH1D("numEta_120","numEta_120",NEtaBins,etaMin,etaMax);
TH1D *numPhi_120 = new TH1D("numPhi_120","numPhi_120",NPhiBins,phiMin,phiMax);
TH1D *r_120 = new TH1D("r_120","r_120",NPtBins,PtMin,PtMax);



// ==========================================================
// main: trigger turn on curves for HLT electrons in MC
// ==========================================================

void triggerAnalysis_PbPb_pho_prompt(
    string inputForest = "/eos/cms/store/group/phys_heavyions/nbarnett/Forests/SpurPbPb2025_RawPrimes/CRAB_UserFiles/crab_forest_SpurPbPb_RPs_11_11_2025_v1/251111_105132/0000/",
    string inputHLT = "", //"/afs/cern.ch/user/k/kdeverea/HLTClayton/CMSSW_15_1_0/src/workstation/HLT_emulation/scripts/PbPb/openHLTfiles/",
    string output_base = "spuriousPbPb_251111",
    int nfiles = 1
  ){

  std::cout << "running triggerAnalysis_ele_mc()" << std::endl;
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
  //HiTree->SetBranchStatus("weight",1);
  HiTree->SetBranchStatus("vz",1);

  HiTree->SetBranchAddress("evt", &evt);
  HiTree->SetBranchAddress("lumi", &lumi);
  HiTree->SetBranchAddress("run", &run);
  HiTree->SetBranchAddress("hiBin", &hiBin);
  HiTree->SetBranchAddress("hiHF", &hiHF);
  //HiTree->SetBranchAddress("weight", &pthat_weight);
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

  HltTree->Print();
  HltTree->Show(0);

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
  //EventTree->SetBranchStatus("pho_genMatchedIndex", 1);
  //EventTree->SetBranchStatus("mcEt", 1);
  //EventTree->SetBranchStatus("mcPID",1);
  //EventTree->SetBranchStatus("mcEta",1);
  //EventTree->SetBranchStatus("mcPhi",1);

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

  //EventTree->SetBranchAddress("pho_genMatchedIndex", &pho_genMatchedIndex);
  //EventTree->SetBranchAddress("mcEt", &mcEt);
  //EventTree->SetBranchAddress("mcPID", &mcPID);
  //EventTree->SetBranchAddress("mcEta", &mcEta);
  //EventTree->SetBranchAddress("mcPhi", &mcPhi);

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


  //ggHiNtuplizer ggHiNtuple;
  //ggHiNtuple.setupTreeForReading(EventTree);

  std::cout << "done" << std::endl;
  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop  =================

  // loop through reco objects
  int n_reco = HltTree-> GetEntries();
  int nphos = 0;
  int npass_leading = 0;
  int npass_trkrcut = 0;
  int npass_l1 = 0;
  int npass_trigger = 0;
  for (ULong64_t i_event = 0; i_event < HltTree->GetEntries(); ++i_event){

    if(i_event%(HltTree->GetEntries()/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    HiTree->GetEntry(i_event);
    HltTree->GetEntry(i_event);
    EventTree->GetEntry(i_event);
  	
    // event cuts
    if(fabs(vz)>15.0) continue;
    //if(hiBin>180) continue;	

    nphos += nPho;


    // ========== skip duplicate reco events ===========
    // TODO


    // ============= find leading electron ==============
    // take highest reco pT as leading electron
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


    // ============= apply track rejection =============
    float sumIso = pfpIso3subUE->at(i_leading) + pfcIso3subUE->at(i_leading) + pfnIso3subUE->at(i_leading);

    // barrel cut
    /*
    if (abs(phoSCEta->at(i_leading)) < 1.442) {
      if (sumIso > 11.697505)                              continue;
      if (phoHoverE->at(i_leading) > 0.247665)             continue;
      if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.012186) continue;
      if (pho_swissCrx->at(i_leading) > 0.9)               continue;
      if (fabs(pho_seedTime->at(i_leading)) > 3.0)         continue;
    }
    else continue;
    */

    if (abs(phoSCEta->at(i_leading)) > 1.442) continue;

    // endcap cut
    /*
    else if (abs(phoSCEta->at(i_leading)) > 1.57 && abs(phoSCEta->at(i_leading)) < 2.1) {
      if (sumIso > 20.911811)                              continue;
      if (phoHoverE->at(i_leading) > 0.398866)             continue;
      if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.044998) continue;
      if (pho_swissCrx->at(i_leading) < 0.9)               continue;
      if (fabs(pho_seedTime->at(i_leading)) > 3.0)         continue;
    }
    else continue;
    */
    npass_trkrcut++;

    // =============== fill histograms ================

    float weight = pthat_weight;

    // fill denominator histograms
    denom->Fill(maxPt, weight);
    
    if(L1_SingleEG7) {
      l1_7->Fill(maxPt, weight);
      npass_l1 ++;
    }

    if(L1_SingleEG15) l1_15->Fill(maxPt, weight);
    if(L1_SingleEG21) l1_21->Fill(maxPt, weight);

    // fill numerator histograms

    // match leading reco electron to HLT emulated object
    /*
    int minDR_HLT = 9999.0;
    for (Int_t i_track = 0; i_track < (int)Ele20Gsf_pt->size(); i_track++){

      float dEta_HLT = eleEta->at(i_leading) - Ele20Gsf_eta->at(i_track);
      float dPhi_HLT = elePhi->at(i_leading) - Ele20Gsf_phi->at(i_track);
      if (dPhi_HLT > TMath::Pi()) dPhi_HLT -= 2*TMath::Pi();
      float dR_HLT = sqrt( dEta_HLT*dEta_HLT + dPhi_HLT*dPhi_HLT );

      if (dR_HLT < minDR_HLT) {
        minDR_HLT = dR_HLT;
      }
    }
    //cout<< "minDR_HLT = " << minDR_HLT << endl;
    if(minDR_HLT > 0.5) continue;
    */
    //npass_hltmatch++;

    if(HLT_HIGEDPhoton20) {
      num_20->Fill(maxPt, weight);
      npass_trigger++;
    }

    if(HLT_HIGEDPhoton30) num_30->Fill(maxPt, weight);

    if(HLT_HIGEDPhoton40) num_40->Fill(maxPt, weight);

    if(HLT_HIGEDPhoton50) num_50->Fill(maxPt, weight);

  } // end of event loop

  std::cout << std::endl;
  std::cout << "nphos = " << nphos  << std::endl;
  std::cout << "npass_leading = " << npass_leading  << std::endl;
  std::cout << "npass_trkrcut = " << npass_trkrcut  << std::endl;
  std::cout << "npass_l1 = " << npass_l1  << std::endl;
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

  r_20->SetLineColor(kRed-4);
  r_30->SetLineColor(kBlue-4);
  r_40->SetLineColor(kGreen+2);
  r_50->SetLineColor(kMagenta-9);
  r_120->SetLineColor(kPink+6);

  r_20->SetMarkerColor(kRed-4);
  r_30->SetMarkerColor(kBlue-4);
  r_40->SetMarkerColor(kGreen+2);
  r_50->SetMarkerColor(kMagenta-9);
  r_120->SetMarkerColor(kPink+6);

  double line_width = 1.8;
  r_20->SetLineWidth(line_width);
  r_30->SetLineWidth(line_width);
  r_40->SetLineWidth(line_width);
  r_50->SetLineWidth(line_width);
  r_120->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_20->SetMarkerSize(marker_size);
  r_30->SetMarkerSize(marker_size);
  r_40->SetMarkerSize(marker_size);
  r_50->SetMarkerSize(marker_size);
  r_120->SetMarkerSize(marker_size);
    
  r_20->SetMarkerStyle(20);
  r_30->SetMarkerStyle(21);
  r_40->SetMarkerStyle(22);
  r_50->SetMarkerStyle(23);
  r_120->SetMarkerStyle(34);

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

  // =============== Fig 1 HLT (no L1 emulation) ==================

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
  //leg->AddEntry(r_50,"HLT_PFJet110_v16");
  //leg->AddEntry(r_120,"HLT_PFJet140_v35"); 
  //leg->SetBorderSize(0);
  r_20->Draw();
  leg->Draw();
  r_30->Draw("same");
  r_40->Draw("same");
  r_50->Draw("same");
  //r_120->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet QCDPhoton");
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  //la->DrawLatexNDC(0.6,0.69,"2025 Run 3 MC");
  la->DrawLatexNDC(0.6,0.69,"|#eta| < 1.442");

  c1->SaveAs(Form("plots/TriggerEfficiency_pho_%s.png", output_base.c_str()));


  // =============== Fig 2 L1  ==================

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
  la2->DrawLatexNDC(0.6,0.69,"|#eta| < 1.442");

  c2->SaveAs(Form("plots/TriggerEfficiency_L1_pho_%s.png", output_base.c_str()));


  auto wf = TFile::Open(Form("output_pho_%s.root", output_base.c_str()),"recreate");

  denom->Write();

  denomEta_40->Write();
  denomEta_60->Write();
  denomEta_70->Write();
  denomEta_80->Write();
  denomEta_100->Write();
  denomEta_120->Write();

  denomPhi_40->Write();
  denomPhi_60->Write();
  denomPhi_70->Write();
  denomPhi_80->Write();
  denomPhi_100->Write();
  denomPhi_120->Write();

  num_20->Write();
  num_30->Write();
  num_40->Write();
  num_50->Write();
  num_120->Write();

  numEta_40->Write();
  numEta_60->Write();
  numEta_70->Write();
  numEta_100->Write();
  numEta_120->Write();

  numPhi_40->Write();
  numPhi_60->Write();
  numPhi_70->Write();
  numPhi_100->Write();
  numPhi_120->Write();

  l1_7->Write();
  l1_15->Write();
  l1_21->Write();

  r_20->Write();
  r_30->Write();
  r_40->Write();
  r_40->Write();
  r_50->Write();
  r_120->Write();

  rl1_7->Write();
  rl1_15->Write();
  rl1_21->Write();

  wf->Close();

}


