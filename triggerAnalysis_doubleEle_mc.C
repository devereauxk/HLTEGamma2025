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
TH1D *denom_endcap = new TH1D("denom_endcap","denom_endcap",NPtBins,PtMin,PtMax);

TH1D *l1_5 = new TH1D("l1_5","l1_5",NPtBins,PtMin,PtMax);
TH1D *l1_5_subleading = new TH1D("l1_5_subleading","l1_5_subleading",NPtBins,PtMin,PtMax);

TH1D *num_10_10 = new TH1D("num_10_10","num_10_10",NPtBins,PtMin,PtMax);
TH1D *num_15_10 = new TH1D("num_15_10","num_15_10",NPtBins,PtMin,PtMax);
TH1D *num_15_15 = new TH1D("num_15_15","num_15_15",NPtBins,PtMin,PtMax);

TH1D *num_10_10_subleading = new TH1D("num_10_10_subleading","num_10_10_subleading",NPtBins,PtMin,PtMax);
TH1D *num_15_10_subleading = new TH1D("num_15_10_subleading","num_15_10_subleading",NPtBins,PtMin,PtMax);
TH1D *num_15_15_subleading = new TH1D("num_15_15_subleading","num_15_15_subleading",NPtBins,PtMin,PtMax);

TH1D *r_10_10 = new TH1D("r_10_10","r_10_10",NPtBins,PtMin,PtMax);
TH1D *r_15_10 = new TH1D("r_15_10","r_15_10",NPtBins,PtMin,PtMax);
TH1D *r_15_15 = new TH1D("r_15_15","r_15_15",NPtBins,PtMin,PtMax);

TH1D *r_10_10_subleading = new TH1D("r_10_10_subleading","r_10_10_subleading",NPtBins,PtMin,PtMax);
TH1D *r_15_10_subleading = new TH1D("r_15_10_subleading","r_15_10_subleading",NPtBins,PtMin,PtMax);
TH1D *r_15_15_subleading = new TH1D("r_15_15_subleading","r_15_15_subleading",NPtBins,PtMin,PtMax);

TH1D *h_hiBin = new TH1D("h_hiBin","h_hiBin",250,0,250);

TH2D *h2_missedEtaPhi = new TH2D("h2_missedEtaPhi","h2_missedEtaPhi",NPhiBins,phiMin,phiMax,NEtaBins,etaMin,etaMax);



// ==========================================================
// main: trigger turn on curves for HLT double electrons in MC
// ==========================================================

void triggerAnalysis_doubleEle_mc(

    // ZtoEE sample
    //string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_Zee_quick/251023_023035/0000/",
    //string inputHLT = "", //"/afs/cern.ch/user/k/kdeverea/HLTClayton/CMSSW_15_1_0/src/workstation/HLT_emulation/scripts/PbPb/openHLTfiles/",
    //string output_base = "MCZee0_100",

    // JpsiToEE sample
    string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000/",
    string inputHLT = "",
    string output_base = "MCJpsiToEE0_100",

    int nfiles = -1,
    float minHiBin = 0.0,
    float maxHiBin = 200.0
  ){

  std::cout << "running triggerAnalysis_doubleEle_mc()" << std::endl;
  std::cout << "input forest directory = " << inputForest  << std::endl;
  std::cout << "input HLT directory    = " << inputHLT  << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  TChain *Tree_HLT_HIEle20Gsf = new TChain("hltobject/HLT_HIEle20Gsf_v");
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
    Tree_HLT_HIEle20Gsf     ->Add(this_hlt_name.c_str());
    HltTree_emulated        ->Add(this_hlt_name.c_str());
  }
  else {
    for(int i = 0; i < nfiles; i++) {
      string this_hlt_name = inputHLT + "/20251030_openHLT_Zee_" + to_string(i) + ".root";
      Tree_HLT_HIEle20Gsf   ->Add(this_hlt_name.c_str());
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
  Int_t           L1_DoubleEG5_BptxAND;
  Int_t           HLT_HIDoubleEle10Gsf;
  Int_t           HLT_HIEle15Ele10Gsf;
  Int_t           HLT_HIDoubleEle15Gsf;
  ULong64_t       forest_Event;
  Int_t           forest_LumiBlock, forest_Run;

  HltTree->SetBranchStatus("*",0);     // disable all branches
  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  HltTree->SetBranchStatus("L1_DoubleEG5_BptxAND", 1);

  HltTree->SetBranchStatus("HLT_HIDoubleEle10Gsf_v16", 1);
  HltTree->SetBranchStatus("HLT_HIEle15Ele10Gsf_v16", 1);
  HltTree->SetBranchStatus("HLT_HIDoubleEle15Gsf_v16", 1);

  HltTree->SetBranchAddress("Event", &forest_Event);
  HltTree->SetBranchAddress("LumiBlock", &forest_LumiBlock);
  HltTree->SetBranchAddress("Run", &forest_Run);

  HltTree->SetBranchAddress("L1_DoubleEG5_BptxAND", &L1_DoubleEG5_BptxAND);

  HltTree->SetBranchAddress("HLT_HIDoubleEle10Gsf_v16", &HLT_HIDoubleEle10Gsf);
  HltTree->SetBranchAddress("HLT_HIEle15Ele10Gsf_v16", &HLT_HIEle15Ele10Gsf);
  HltTree->SetBranchAddress("HLT_HIDoubleEle15Gsf_v16", &HLT_HIDoubleEle15Gsf);

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

  // HLT_HIEle20Gsf from emulated HLT
  vector<float>           *Ele20Gsf_pt = nullptr;
  vector<float>           *Ele20Gsf_eta = nullptr;
  vector<float>           *Ele20Gsf_phi = nullptr;

  Tree_HLT_HIEle20Gsf->SetBranchStatus("*",0);     // disable all branches
  Tree_HLT_HIEle20Gsf->SetBranchStatus("pt", 1);
  Tree_HLT_HIEle20Gsf->SetBranchStatus("eta", 1);
  Tree_HLT_HIEle20Gsf->SetBranchStatus("phi", 1);

  Tree_HLT_HIEle20Gsf->SetBranchAddress("pt", &Ele20Gsf_pt);
  Tree_HLT_HIEle20Gsf->SetBranchAddress("eta", &Ele20Gsf_eta);
  Tree_HLT_HIEle20Gsf->SetBranchAddress("phi", &Ele20Gsf_phi);
  */


  // ================= Event Tree =================
  vector<float>   *elePt = nullptr;
  vector<float>   *eleEta = nullptr;
  vector<float>   *elePhi = nullptr;
  Int_t           nEle;
  vector<float>   *eleHoverE = nullptr;
  vector<float>   *eleSigmaIEtaIEta_2012 = nullptr;
  vector<float>   *eledEtaSeedAtVtx = nullptr;
  vector<float>   *eledPhiAtVtx = nullptr;
  vector<float>   *eleHoverEBc = nullptr;
  vector<float>   *eleEoverPInv = nullptr;
  vector<int>     *eleMissHits = nullptr;
  vector<float>   *eleIP3D = nullptr;
  vector<float>   *eleSCEta = nullptr;
  vector<float>   *eleSCRawEn = nullptr;
  vector<int>     *mcPID = nullptr;
  vector<int>     *mcMomPID = nullptr;
  vector<float>   *mcPt = nullptr;
  vector<float>   *mcEta = nullptr;
  vector<float>   *mcPhi = nullptr;

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

  // mc branches
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

  EventTree->SetBranchStatus("eleSigmaIEtaIEta_2012",1);
  EventTree->SetBranchStatus("eledEtaSeedAtVtx",1);
  EventTree->SetBranchStatus("eledPhiAtVtx",1);
  EventTree->SetBranchStatus("eleHoverEBc",1);
  EventTree->SetBranchStatus("eleEoverPInv",1);
  EventTree->SetBranchStatus("eleMissHits",1);
  EventTree->SetBranchStatus("eleIP3D",1);
  EventTree->SetBranchStatus("eleSCEta",1);
  EventTree->SetBranchStatus("eleSCRawEn",1);

  EventTree->SetBranchAddress("mcPID", &mcPID);
  EventTree->SetBranchAddress("mcMomPID", &mcMomPID);
  EventTree->SetBranchAddress("mcPt", &mcPt);
  EventTree->SetBranchAddress("mcEta", &mcEta);
  EventTree->SetBranchAddress("mcPhi", &mcPhi);

  EventTree->SetBranchAddress("nEle", &nEle);
  EventTree->SetBranchAddress("elePt", &elePt);
  EventTree->SetBranchAddress("eleEta", &eleEta);
  EventTree->SetBranchAddress("eleHoverE", &eleHoverE);
  EventTree->SetBranchAddress("elePhi", &elePhi);

  EventTree->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
  EventTree->SetBranchAddress("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
  EventTree->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
  EventTree->SetBranchAddress("eleHoverEBc", &eleHoverEBc);
  EventTree->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
  EventTree->SetBranchAddress("eleMissHits", &eleMissHits);
  EventTree->SetBranchAddress("eleIP3D", &eleIP3D);
  EventTree->SetBranchAddress("eleSCEta", &eleSCEta);
  EventTree->SetBranchAddress("eleSCRawEn", &eleSCRawEn);
  

  //ggHiNtuplizer ggHiNtuple;
  //ggHiNtuple.setupTreeForReading(EventTree);

  std::cout << "done" << std::endl;
  //Long64_t entriesHLT = HltTree_emulated->GetEntries();
  //std::cout << "HLT entries = " << entriesHLT << std::endl;

  // fill event matcher with HLT events
  /*
  EventMatcher* emTrig = new EventMatcher();
  for (ULong64_t j_event = 0; j_event < HltTree_emulated->GetEntries(); ++j_event){ 

    HltTree_emulated->GetEntry(j_event);

    emTrig->addEvent(hlt_Run, hlt_LumiBlock, hlt_Event, j_event);
  }
  */

  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop  =================

  // loop through reco objects
  int n_reco = HltTree-> GetEntries();
  int neles = 0;
  int npass_emmatch = 0;
  int npass_leading = 0;
  int npass_trkrcut = 0;
  int npass_genmatch = 0;
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

    neles += nEle;

    // ========== skip duplicate reco events ===========
    // TODO

    // ===== find HLT emulation match of reco event =====
    // TODO

    npass_emmatch++;

    // ============= find leading electron ==============
    // take highest reco pT as leading electron
    float maxPt = 0;
    int i_leading = -1;
    float maxPt_subleading = 0;
    int i_subleading = -1;
    
    // loop over reco
    for(Int_t i_track = 0; i_track < nEle; i_track++){

      // compare with previous leading candidate
      if(elePt->at(i_track) > maxPt) {
        maxPt = elePt->at(i_track);
        i_leading = i_track;
      } else if(elePt->at(i_track) > maxPt_subleading) {
        maxPt_subleading = elePt->at(i_track);
        i_subleading = i_track;
      }

    }
    if (i_leading == -1 || i_subleading == -1) continue;
    npass_leading++;


    // ============= apply track rejection =============

    bool isBarrel_leading = false;
    bool isEndcap_leading = false;
    bool isBarrel_subleading = false;
    bool isEndcap_subleading = false;

    // barrel cut
    if (abs(eleEta->at(i_leading)) < 1.442) {
      if (eleSigmaIEtaIEta_2012->at(i_leading) > 0.012)   continue;
      if (eledEtaSeedAtVtx->at(i_leading) > 0.0037)       continue;
      if (eledPhiAtVtx->at(i_leading) > 0.1280)           continue;
      if (eleEoverPInv->at(i_leading) > 0.1065)           continue;
      isBarrel_leading = true;
    }

    // endcap cut
    else if (abs(eleEta->at(i_leading)) > 1.556 && abs(eleEta->at(i_leading)) < 2.1) {
      if (eleSigmaIEtaIEta_2012->at(i_leading) > 0.0376)  continue;
      if (eledEtaSeedAtVtx->at(i_leading) > 0.0074)       continue;
      if (eledPhiAtVtx->at(i_leading) > 0.2085)           continue;
      if (eleEoverPInv->at(i_leading) > 0.1138)           continue;
      isEndcap_leading = true;
    }
    else continue;
    if ( eleMissHits->at(i_leading) > 3 )         continue;
    if ( eleIP3D->at(i_leading) >= 0.03 )         continue;

    // barrel cut
    if (abs(eleEta->at(i_subleading)) < 1.442) {
      if (eleSigmaIEtaIEta_2012->at(i_subleading) > 0.012)   continue;
      if (eledEtaSeedAtVtx->at(i_subleading) > 0.0037)       continue;
      if (eledPhiAtVtx->at(i_subleading) > 0.1280)           continue;
      if (eleEoverPInv->at(i_subleading) > 0.1065)           continue;
      isBarrel_subleading = true;
    }

    // endcap cut
    else if (abs(eleEta->at(i_subleading)) > 1.556 && abs(eleEta->at(i_subleading)) < 2.1) {
      if (eleSigmaIEtaIEta_2012->at(i_subleading) > 0.0376)  continue;
      if (eledEtaSeedAtVtx->at(i_subleading) > 0.0074)       continue;
      if (eledPhiAtVtx->at(i_subleading) > 0.2085)           continue;
      if (eleEoverPInv->at(i_subleading) > 0.1138)           continue;
      isEndcap_subleading = true;
    }
    else continue;
    if ( eleMissHits->at(i_subleading) > 3 )         continue;
    if ( eleIP3D->at(i_subleading) >= 0.03 )         continue;
    

    npass_trkrcut++;


    // =============== match to gen level ===============
    float dRMax = 0.0025;
    float relativeDiffpTMax = 0.1;

    int i_genleading = -1;
    float maxPt_gen = 0;
    int i_genmatch = -1;
    float maxPt_genmatch = 0;

    int i_gensubleading = -1;
    float maxPt_gensubleading = 0;
    int i_gensubmatch = -1;
    float maxPt_gensubmatch = 0;
    
    // loop over gen
    // 1) find leading gen-level electron
    // 2) find best gen-level electron match to leading reco-level electron
    for(Int_t j_track = 0; j_track < (int)mcPID->size(); j_track++){

      // gen electron selection
      if( abs(mcPID->at(j_track)) != 11 ) continue;

      // compare with previous leading candidate
      if( mcPt->at(j_track) > maxPt_gen ) {
        maxPt_gen = mcPt->at(j_track);
        i_genleading = j_track;
      } else if( mcPt->at(j_track) > maxPt_gensubleading ) {
        maxPt_gensubleading = mcPt->at(j_track);
        i_gensubleading = j_track;
      }

      // gen-reco matching
      // deltaR requirement
      float dEta = eleEta->at(i_leading) - mcEta->at(j_track);
      float dPhi = elePhi->at(i_leading) - mcPhi->at(j_track);
      float dR = sqrt( dEta*dEta + dPhi*dPhi );

      float relativeDiffpT = fabs( elePt->at(i_leading) - mcPt->at(j_track) ) / mcPt->at(j_track);

      if (dR < dRMax && relativeDiffpT < relativeDiffpTMax) {  // matched candidate
        // compare with previous leading candidate
        if( mcPt->at(j_track) > maxPt_genmatch ) { // find the leading mcEt in the event
          maxPt_genmatch = mcPt->at(j_track);
          i_genmatch = j_track;
        }
      }

      float dEta_sub = eleEta->at(i_subleading) - mcEta->at(j_track);
      float dPhi_sub = elePhi->at(i_subleading) - mcPhi->at(j_track);
      float dR_sub = sqrt( dEta_sub*dEta_sub + dPhi_sub*dPhi_sub );

      float relativeDiffpT_sub = fabs( elePt->at(i_subleading) - mcPt->at(j_track) ) / mcPt->at(j_track);

      if (dR_sub < dRMax && relativeDiffpT_sub < relativeDiffpTMax) {  // matched candidate
        // compare with previous leading candidate
        if( mcPt->at(j_track) > maxPt_gensubmatch ) { // find the leading mcEt in the event
          maxPt_gensubmatch = mcPt->at(j_track);
          i_gensubmatch = j_track;
        }
      }

    }
    if (i_genleading == -1 || i_genleading != i_genmatch) continue;
    if (i_gensubleading == -1 || i_gensubleading != i_gensubmatch) continue;

    npass_genmatch++;


    // =============== fill histograms ================

    float weight = 1;

    // fill denominator histograms
    denom->Fill(maxPt, weight);

    if(L1_DoubleEG5_BptxAND) {
      l1_5->Fill(maxPt, weight);
      l1_5_subleading->Fill(maxPt_subleading, weight);
      npass_l1++;
    }

    // fill numerator histograms
    if(HLT_HIDoubleEle10Gsf) {
      num_10_10->Fill(maxPt, weight);
      num_10_10_subleading->Fill(maxPt_subleading, weight);
      npass_trigger++;
    }

    if(HLT_HIEle15Ele10Gsf) {
      num_15_10->Fill(maxPt, weight);
      num_15_10_subleading->Fill(maxPt_subleading, weight);
    }

    if(HLT_HIDoubleEle15Gsf) {
      num_15_15->Fill(maxPt, weight);
      num_15_15_subleading->Fill(maxPt_subleading, weight);
    }

    //if(L1_SingleEG5 && !HLT_HIEle20Gsf) {
    //  h2_missedEtaPhi->Fill(elePhi->at(i_leading), eleEta->at(i_leading));
    //}


  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIEle20Gsf results ... " << std::endl;
  std::cout << "neles = " << neles << std::endl;
  std::cout << "npass_emmatch = " << npass_emmatch << std::endl;
  std::cout << "npass_leading = " << npass_leading  << std::endl;
  std::cout << "npass_trkrcut = " << npass_trkrcut  << std::endl;
  std::cout << "npass_genmatch = " << npass_genmatch  << std::endl;
  std::cout << "npass_trigger = " << npass_trigger  << std::endl;

  r_10_10->Divide(num_10_10,l1_5,1,1,"B");
  r_15_10->Divide(num_15_10,l1_5,1,1,"B");
  r_15_15->Divide(num_15_15,l1_5,1,1,"B");

  r_10_10_subleading->Divide(num_10_10_subleading,l1_5_subleading,1,1,"B");
  r_15_10_subleading->Divide(num_15_10_subleading,l1_5_subleading,1,1,"B");
  r_15_15_subleading->Divide(num_15_15_subleading,l1_5_subleading,1,1,"B");

  r_10_10->SetLineColor(kRed-4);
  r_15_10->SetLineColor(kBlue-4);
  r_15_15->SetLineColor(kGreen+2);
  r_10_10_subleading->SetLineColor(kRed-4);
  r_15_10_subleading->SetLineColor(kBlue-4);
  r_15_15_subleading->SetLineColor(kGreen+2);

  r_10_10->SetMarkerColor(kRed-4);
  r_15_10->SetMarkerColor(kBlue-4);
  r_15_15->SetMarkerColor(kGreen+2);
  r_10_10_subleading->SetMarkerColor(kRed-4);
  r_15_10_subleading->SetMarkerColor(kBlue-4);
  r_15_15_subleading->SetMarkerColor(kGreen+2);

  double line_width = 1.8;
  r_10_10->SetLineWidth(line_width);
  r_15_10->SetLineWidth(line_width);
  r_15_15->SetLineWidth(line_width);
  r_10_10_subleading->SetLineWidth(line_width);
  r_15_10_subleading->SetLineWidth(line_width);
  r_15_15_subleading->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_10_10->SetMarkerSize(marker_size);
  r_15_10->SetMarkerSize(marker_size);
  r_15_15->SetMarkerSize(marker_size);
  r_10_10_subleading->SetMarkerSize(marker_size);
  r_15_10_subleading->SetMarkerSize(marker_size);
  r_15_15_subleading->SetMarkerSize(marker_size);
    
  r_10_10->SetMarkerStyle(20);
  r_15_10->SetMarkerStyle(21);
  r_15_15->SetMarkerStyle(22);
  r_10_10_subleading->SetMarkerStyle(20);
  r_15_10_subleading->SetMarkerStyle(21);
  r_15_15_subleading->SetMarkerStyle(22);


  // =============== Fig HLT (no L1 emulation) LEADING ELE ==================

  TCanvas *c1 = new TCanvas("c1","c1",700,600);
  c1->cd();
  TPad *p1 = new TPad("p1","p1",0,0,1,1);
  p1->SetLeftMargin(0.13);
  p1->SetBottomMargin(0.14);
  p1->Draw();
  p1->cd();
  r_10_10->SetTitle("");
  r_10_10->SetStats(0);
  r_10_10->GetXaxis()->SetTitleSize(0.05);
  r_10_10->GetYaxis()->SetTitleSize(0.05);
  r_10_10->GetXaxis()->SetTitle("leading electron #font[52]{p}_{T} [GeV]");
  r_10_10->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_10_10,"HLT_HIDoubleEle10Gsf_v16");
  leg->AddEntry(r_15_10,"HLT_HIEle15Ele10Gsf_v16");
  leg->AddEntry(r_15_15,"HLT_HIDoubleEle15Gsf_v16");
  //leg->SetBorderSize(0);
  r_10_10->Draw();
  leg->Draw();
  r_15_10->Draw("same");
  r_15_15->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la->DrawLatexNDC(0.6,0.65,"Barrel + Endcap");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("plots/TriggerEfficiency_doubleEle_%s_leading.png", output_base.c_str()));


  // =============== Fig HLT (no L1 emulation) SUBLEADING ELE ==================
  TCanvas *c1_endcap = new TCanvas("c1_endcap","c1_endcap",700,600);
  c1_endcap->cd();
  TPad *p1_endcap = new TPad("p1_endcap","p1_endcap",0,0,1,1);
  p1_endcap->SetLeftMargin(0.13);
  p1_endcap->SetBottomMargin(0.14);
  p1_endcap->Draw();
  p1_endcap->cd();
  r_10_10_subleading->SetTitle("");
  r_10_10_subleading->SetStats(0);
  r_10_10_subleading->GetXaxis()->SetTitleSize(0.05);
  r_10_10_subleading->GetYaxis()->SetTitleSize(0.05);
  r_10_10_subleading->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_10_10_subleading->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend* leg_subleading = new TLegend(0.55,0.3,0.88,0.5);
  leg_subleading->AddEntry(r_10_10_subleading,"HLT_HIDoubleEle10Gsf_v16");
  leg_subleading->AddEntry(r_15_10_subleading,"HLT_HIEle15Ele10Gsf_v16");
  leg_subleading->AddEntry(r_15_15_subleading,"HLT_HIDoubleEle15Gsf_v16");
  //leg_subleading->SetBorderSize(0);
  r_10_10_subleading->Draw();
  leg_subleading->Draw();
  r_15_10_subleading->Draw("same");
  r_15_15_subleading->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la_endcap->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la_endcap->DrawLatexNDC(0.6,0.65,"Barrel + Endcap");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("plots/TriggerEfficiency_doubleEle_%s_subleading.png", output_base.c_str()));


  // =============== Fig L1 ==================

  /*
  TCanvas *c2 = new TCanvas("c2","c2",700,600);
  c2->cd();
  TPad *p2 = new TPad("p2","p2",0,0,1,1);
  p2->SetLeftMargin(0.13);
  p2->SetBottomMargin(0.14);
  p2->Draw();
  p2->cd();
  rl1_5->SetTitle("");
  rl1_5->SetStats(0);
  rl1_5->GetXaxis()->SetTitleSize(0.05);
  rl1_5->GetYaxis()->SetTitleSize(0.05);
  rl1_5->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  rl1_5->GetYaxis()->SetTitle("L1 Trigger efficiency");
  TLegend *leg2 = new TLegend(0.55,0.3,0.88,0.5);
  leg2->AddEntry(rl1_5,"L1_SingleEG5_BptxAND");
  //leg2->SetBorderSize(0);
  rl1_5->Draw();
  leg2->Draw();
  rl1_15->Draw("same");
  rl1_21->Draw("same");

  TLatex *la2 = new TLatex();
  la2->SetTextFont(42);
  la2->SetTextSize(0.03);

  la2->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la2->DrawLatexNDC(0.72,0.92,"L1 emulation");
  la2->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la2->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c2->SaveAs(Form("plots/TriggerEfficiency_ele_%s_L1.png", output_base.c_str()));


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
  c3->SaveAs(Form("plots/hiBin_ele_%s.png", output_base.c_str()));


  // ========== missed eta-phi plot ==========
  TCanvas *c4 = new TCanvas("c4","c4",700,600);
  c4->cd();
  TPad *p4 = new TPad("p4","p4",0,0,1,1);
  p4->SetLeftMargin(0.13);
  p4->SetBottomMargin(0.14);
  p4->Draw();
  p4->cd();
  h2_missedEtaPhi->SetTitle("");
  h2_missedEtaPhi->SetStats(0);
  h2_missedEtaPhi->GetXaxis()->SetTitleSize(0.05);
  h2_missedEtaPhi->GetYaxis()->SetTitleSize(0.05);
  h2_missedEtaPhi->GetXaxis()->SetTitle("phi");
  h2_missedEtaPhi->GetYaxis()->SetTitle("eta");
  h2_missedEtaPhi->Draw("COLZ");
  c4->SaveAs(Form("plots/MissedEtaPhi_ele_%s.png", output_base.c_str()));
  


  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/output_ele_%s.root", output_base.c_str()),"recreate");

  denom->Write();
  denom_endcap->Write();

  num_10_10->Write();
  num_15_10->Write();
  num_15_15->Write();
  num_50->Write();

  l1_5->Write();
  l1_15->Write();
  l1_21->Write();

  num_10_10_endcap->Write();
  num_15_10_endcap->Write();
  num_15_15_endcap->Write();
  num_50_endcap->Write();

  l1_5_endcap->Write();
  l1_15_endcap->Write();
  l1_21_endcap->Write();

  r_10_10->Write();
  r_30->Write();
  r_40->Write();
  r_50->Write();

  r_10_10_subleading->Write();
  r_30_endcap->Write();
  r_40_endcap->Write();
  r_50_endcap->Write();

  wf->Close();
  */

}


