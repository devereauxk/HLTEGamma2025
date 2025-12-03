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

TH2D *h2_missedEtaPhi = new TH2D("h2_missedEtaPhi","h2_missedEtaPhi",NPhiBins,phiMin,phiMax,NEtaBins,etaMin,etaMax);


// debugging cuts
map<string, TH1D*> h_cuts_pass = {
  {"eleMissHits",           new TH1D("eleMissHits","eleMissHits",10,0,4)},
  {"eleIP3D",               new TH1D("eleIP3D","eleIP3D",10,0,0.03)},
  {"eleSigmaIEtaIEta_2012", new TH1D("eleSigmaIEtaIEta_2012","eleSigmaIEtaIEta_2012",10,0,0.02)},
  {"eledEtaSeedAtVtx",      new TH1D("eledEtaSeedAtVtx","eledEtaSeedAtVtx",20,-0.004,0.004)},
  {"eledPhiAtVtx",          new TH1D("eledPhiAtVtx","eledPhiAtVtx",10,0,0.08)},
  {"eleEoverPInv",          new TH1D("eleEoverPInv","eleEoverPInv",10,0,0.04)},
  {"eleHoverE",             new TH1D("eleHoverE","eleHoverE",10,0,0.14)},
  {"deltaR",                new TH1D("deltaR","deltaR",10,0,0.01)}
};

map<string, TH1D*> h_cuts_fail = {
  {"eleMissHits",           new TH1D("eleMissHits_fail","eleMissHits_fail",10,0,4)},
  {"eleIP3D",               new TH1D("eleIP3D_fail","eleIP3D_fail",10,0,0.03)},
  {"eleSigmaIEtaIEta_2012", new TH1D("eleSigmaIEtaIEta_2012_fail","eleSigmaIEtaIEta_2012_fail",10,0,0.02)},
  {"eledEtaSeedAtVtx",      new TH1D("eledEtaSeedAtVtx_fail","eledEtaSeedAtVtx_fail",20,-0.004,0.004)},
  {"eledPhiAtVtx",          new TH1D("eledPhiAtVtx_fail","eledPhiAtVtx_fail",10,0,0.08)},
  {"eleEoverPInv",          new TH1D("eleEoverPInv_fail","eleEoverPInv_fail",10,0,0.04)},
  {"eleHoverE",             new TH1D("eleHoverE_fail","eleHoverE_fail",10,0,0.14)},
  {"deltaR",                new TH1D("deltaR_fail","deltaR_fail",10,0,0.01)}
};


unsigned long long keyFromRunLumiEvent(Int_t run,
                                       Int_t lumi,
                                       ULong64_t event)
{
  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUNLUMIEVENTKEY WARNING : \'" << event
              << "\' is greated that event limit 10^10. returning key 0"
              << std::endl;
    return key;
  }

  key += runMult* static_cast<unsigned long long>(run);
  key += lumiMult*static_cast<unsigned long long>(lumi);
  key += evtMult*event;

  //std::cout << "key = " << key << std::endl;
  
  return key;
  
}

void plotCutPassAndFail(TH1D *h_pass, TH1D *h_fail, string name) {
  TCanvas *c3 = new TCanvas("c3","c3",700,600);
  c3->cd();
  TPad *p3 = new TPad("p3","p3",0,0,1,1);
  p3->SetLeftMargin(0.13);
  p3->SetBottomMargin(0.14);
  p3->Draw();
  p3->cd();
  h_pass->SetTitle(name.c_str());
  h_pass->SetStats(0);
  h_pass->GetXaxis()->SetTitleSize(0.05);
  h_pass->GetYaxis()->SetTitleSize(0.05);
  h_pass->GetXaxis()->SetTitle(Form("%s value", name.c_str()));
  h_pass->GetYaxis()->SetTitle("Counts");

  TLegend *leg = new TLegend(0.64,0.74,0.84,0.88);
  leg->AddEntry(h_pass,"pass","l");
  leg->AddEntry(h_fail,"fail","l");

  h_pass->Scale(1.0/h_pass->Integral());
  h_fail->Scale(1.0/h_fail->Integral());

  h_pass->SetLineColor(kBlue);
  h_pass->Draw();
  h_fail->SetLineColor(kRed);
  h_fail->Draw("SAME");
  leg->Draw();

  c3->SaveAs(Form("plots/debug_cut_%s.png", name.c_str()));
}


// ==========================================================
// main: trigger turn on curves for HLT electrons in MC
// ==========================================================

void debug_triggerAnalysis_ele_mc(

    // ZtoEE sample
    string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_Zee_quick/251023_023035/0000/",
    string inputHLT = "", //"/afs/cern.ch/user/k/kdeverea/HLTClayton/CMSSW_15_1_0/src/workstation/HLT_emulation/scripts/PbPb/openHLTfiles/",
    string output_base = "MCZee0_100_v2",

    // JpsiToEE sample
    //string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000/",
    //string inputHLT = "",
    //string output_base = "MCJpsiToEE",

    int nfiles = 200,
    float minHiBin = 0.0,
    float maxHiBin = 200.0
  ){

  std::cout << "running triggerAnalysis_ele_mc()" << std::endl;
  std::cout << "input forest directory = " << inputForest  << std::endl;
  std::cout << "input HLT directory    = " << inputHLT  << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  //TChain *Tree_HLT_HIEle20Gsf = new TChain("hltobject/HLT_HIEle20Gsf_v");
  //TChain *HltTree_emulated = new TChain("hltanalysis/HltTree");
  
  std::cout<< "Adding input files...";

  // forest files
  if (nfiles == -1) {
    string this_forest_name = inputForest + "/merged.root";
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

  TFile *f_correction = new TFile("output/corrections.root", "READ");
  TH1D* correction = (TH1D*) f_correction->Get("correction_eledEtaSeedAtVtx");

  // hlt files
  //string this_hlt_name = inputHLT + "/*openHLT_Zee_*.root";
  //Tree_HLT_HIEle20Gsf     ->Add(this_hlt_name.c_str());
  //HltTree_emulated        ->Add(this_hlt_name.c_str());

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

  HltTree->SetBranchStatus("HLT_HIEle20Gsf_v16", 1);
  HltTree->SetBranchStatus("HLT_HIEle30Gsf_v16", 1);
  HltTree->SetBranchStatus("HLT_HIEle40Gsf_v16", 1);
  HltTree->SetBranchStatus("HLT_HIEle50Gsf_v16", 1);

  HltTree->SetBranchAddress("Event", &forest_Event);
  HltTree->SetBranchAddress("LumiBlock", &forest_LumiBlock);
  HltTree->SetBranchAddress("Run", &forest_Run);

  HltTree->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7);
  HltTree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15);
  HltTree->SetBranchAddress("L1_SingleEG21_BptxAND", &L1_SingleEG21);

  HltTree->SetBranchAddress("HLT_HIEle20Gsf_v16", &HLT_HIEle20Gsf);
  HltTree->SetBranchAddress("HLT_HIEle30Gsf_v16", &HLT_HIEle30Gsf);
  HltTree->SetBranchAddress("HLT_HIEle40Gsf_v16", &HLT_HIEle40Gsf);
  HltTree->SetBranchAddress("HLT_HIEle50Gsf_v16", &HLT_HIEle50Gsf);

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
  EventTree->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
  EventTree->SetBranchAddress("eleMissHits", &eleMissHits);
  EventTree->SetBranchAddress("eleIP3D", &eleIP3D);
  EventTree->SetBranchAddress("eleSCEta", &eleSCEta);
  EventTree->SetBranchAddress("eleSCRawEn", &eleSCRawEn);
  

  std::cout << "done" << std::endl;
  //Long64_t entriesHLT = HltTree_emulated->GetEntries();
  //std::cout << "HLT entries = " << entriesHLT << std::endl;

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
  int nfail_trigger = 0;
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
    
    // loop over reco
    for(Int_t i_track = 0; i_track < nEle; i_track++){

      // compare with previous leading candidate
      if(elePt->at(i_track) > maxPt) { // find the leading elePt in the event
        maxPt = elePt->at(i_track);
        i_leading = i_track;
      }

    }
    if (i_leading == -1) continue;
    npass_leading++;
    

    // ============= apply track rejection =============
    bool isBarrel = false;
    bool isEndcap = false;

    // barrel cut
    if (abs(eleEta->at(i_leading)) < 1.442) {
      if (eleMissHits->at(i_leading) > 3)                 continue;
      if (eleIP3D->at(i_leading) > 0.03)                  continue;
      if (eleSigmaIEtaIEta_2012->at(i_leading) > 0.012)   continue;
      if (eledEtaSeedAtVtx->at(i_leading) > 0.0037)       continue;
      if (eledPhiAtVtx->at(i_leading) > 0.1280)           continue;
      if (eleEoverPInv->at(i_leading) > 0.1065)           continue;
      if (eleHoverE->at(i_leading) > 0.13)                continue;
      isBarrel = true;
    }

    // endcap cut
    else if (abs(eleEta->at(i_leading)) > 1.556 && abs(eleEta->at(i_leading)) < 2.1) {
      if ( eleMissHits->at(i_leading) > 3 )               continue;
      if ( eleIP3D->at(i_leading) >= 0.03 )               continue;
      if (eleSigmaIEtaIEta_2012->at(i_leading) > 0.0376)  continue;
      if (eledEtaSeedAtVtx->at(i_leading) > 0.0074)       continue;
      if (eledPhiAtVtx->at(i_leading) > 0.2085)           continue;
      if (eleEoverPInv->at(i_leading) > 0.1138)           continue;
      if (eleHoverE->at(i_leading) > 0.13)                continue;
      isEndcap = true;
    }
    else continue;

    npass_trkrcut++;


    // =============== match to gen level ===============
    int i_genleading = -1;
    float maxPt_gen = 0;
    int i_genmatch = -1;
    float maxPt_genmatch = 0;
    float deltaR_genmatch = -1;
    
    // loop over gen
    // 1) find leading gen-level electron
    // 2) find best gen-level electron match to leading reco-level electron
    for(Int_t j_track = 0; j_track < (int)mcPID->size(); j_track++){

      // gen electron selection
      if( abs(mcPID->at(j_track)) != 11) continue;

      // compare with previous leading candidate
      if( mcPt->at(j_track) > maxPt_gen ) {
        maxPt_gen = mcPt->at(j_track);
        i_genleading = j_track;
      }

      // gen-reco matching
      // deltaR requirement
      float dEta = eleEta->at(i_leading) - mcEta->at(j_track);
      float dPhi = elePhi->at(i_leading) - mcPhi->at(j_track);
      float dR = sqrt( dEta*dEta + dPhi*dPhi );

      float relativeDiffpT = fabs( elePt->at(i_leading) - mcPt->at(j_track) ) / mcPt->at(j_track);

      if (dR < 0.0225 && relativeDiffpT < 0.5) {  // matched candidate
        // compare with previous leading candidate
        if( mcPt->at(j_track) > maxPt_genmatch ) { // find the leading mcEt in the event
          maxPt_genmatch = mcPt->at(j_track);
          i_genmatch = j_track;
          deltaR_genmatch = dR;
        }
      }
    }
    if (i_genleading == -1 || i_genleading != i_genmatch) continue;
    //cout<<"i_genleading: "<<i_genleading<<" i_genmatch: "<<i_genmatch<<endl;

    npass_genmatch++;


    // =============== fill histograms ================

    float weight = 1;
    if (!HLT_HIEle30Gsf) {
      weight *= correction->GetBinContent(correction->FindBin(eledEtaSeedAtVtx->at(i_leading)));
      //cout<<"eledEtaSeedAtVtx: "<<eledEtaSeedAtVtx->at(i_leading)<<" weight "<<weight<<endl;
    }

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
      if(HLT_HIEle20Gsf) {
        num_20->Fill(maxPt, weight);
        npass_trigger++;
      }

      if(HLT_HIEle30Gsf) { 
        num_30->Fill(maxPt, weight);

        // fill cut histograms
        if (maxPt > 35) {
          for(auto const& [name, hist] : h_cuts_pass) {
            if(name == "eleMissHits")           hist->Fill(eleMissHits->at(i_leading));
            if(name == "eleIP3D")               hist->Fill(eleIP3D->at(i_leading));
            if(name == "eleSigmaIEtaIEta_2012") hist->Fill(eleSigmaIEtaIEta_2012->at(i_leading));
            if(name == "eledEtaSeedAtVtx")      hist->Fill(eledEtaSeedAtVtx->at(i_leading), weight);
            if(name == "eledPhiAtVtx")          hist->Fill(eledPhiAtVtx->at(i_leading));
            if(name == "eleEoverPInv")          hist->Fill(eleEoverPInv->at(i_leading));
            if(name == "eleHoverE")             hist->Fill(eleHoverE->at(i_leading));
            if(name == "deltaR")                hist->Fill(deltaR_genmatch);
          }
        }
      } else {
        // fill cut histograms
        if (maxPt > 35) {
          for(auto const& [name, hist] : h_cuts_fail) {
            if(name == "eleMissHits")           hist->Fill(eleMissHits->at(i_leading));
            if(name == "eleIP3D")               hist->Fill(eleIP3D->at(i_leading));
            if(name == "eleSigmaIEtaIEta_2012") hist->Fill(eleSigmaIEtaIEta_2012->at(i_leading));
            if(name == "eledEtaSeedAtVtx")      hist->Fill(eledEtaSeedAtVtx->at(i_leading), weight);
            if(name == "eledPhiAtVtx")          hist->Fill(eledPhiAtVtx->at(i_leading));
            if(name == "eleEoverPInv")          hist->Fill(eleEoverPInv->at(i_leading));
            if(name == "eleHoverE")             hist->Fill(eleHoverE->at(i_leading));
            if(name == "deltaR")                hist->Fill(deltaR_genmatch);
          }
        }
        nfail_trigger++;
      }

      if(HLT_HIEle40Gsf) num_40->Fill(maxPt, weight);

      if(HLT_HIEle50Gsf) num_50->Fill(maxPt, weight);

      if(L1_SingleEG7 && !HLT_HIEle20Gsf) {
        h2_missedEtaPhi->Fill(elePhi->at(i_leading), eleEta->at(i_leading));
      }

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
      if(HLT_HIEle20Gsf) {
        num_20_endcap->Fill(maxPt, weight);
      }

      if(HLT_HIEle30Gsf) num_30_endcap->Fill(maxPt, weight);

      if(HLT_HIEle40Gsf) num_40_endcap->Fill(maxPt, weight);

      if(HLT_HIEle50Gsf) num_50_endcap->Fill(maxPt, weight);

    }

  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIEle20Gsf results ... " << std::endl;
  std::cout << "neles = " << neles << std::endl;
  std::cout << "npass_leading = " << npass_leading  << std::endl;
  std::cout << "npass_trkrcut = " << npass_trkrcut  << std::endl;
  std::cout << "npass_genmatch = " << npass_genmatch  << std::endl;
  std::cout << "npass_l1       = " << npass_l1  << std::endl;
  std::cout << "npass_trigger = " << npass_trigger  << std::endl;
  std::cout << "nfail_trigger = " << nfail_trigger  << std::endl;


  // cut debug plots
  TFile *out_correction = new TFile("output/corrections.root", "RECREATE");
  vector<TH1D*> corrections;
  for(auto const& [name, hist] : h_cuts_pass) {
    // plot
    plotCutPassAndFail(hist, h_cuts_fail[name], name);

    // calculate correction, c = pass / fail
    TH1D * this_correction = (TH1D*)hist->Clone(Form("correction_%s", name.c_str()));
    this_correction->Divide(h_cuts_fail[name]);
    this_correction->Write();
  }
  out_correction->Close();

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
  r_20->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_20->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_20,"HLT_HIEle20Gsf_v16");
  leg->AddEntry(r_30,"HLT_HIEle30Gsf_v16");
  leg->AddEntry(r_40,"HLT_HIEle40Gsf_v16");
  leg->AddEntry(r_50,"HLT_HIEle50Gsf_v16");
  //leg->SetBorderSize(0);
  r_20->Draw();
  leg->Draw();
  r_30->Draw("same");
  r_40->Draw("same");
  r_50->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("plots/debug_TriggerEfficiency_ele_%s_barrel.png", output_base.c_str()));


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
  r_20_endcap->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_20_endcap->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg_endcap = new TLegend(0.55,0.3,0.88,0.5);
  leg_endcap->AddEntry(r_20_endcap,"HLT_HIEle20Gsf_v16");
  leg_endcap->AddEntry(r_30_endcap,"HLT_HIEle30Gsf_v16");
  leg_endcap->AddEntry(r_40_endcap,"HLT_HIEle40Gsf_v16");
  leg_endcap->AddEntry(r_50_endcap,"HLT_HIEle50Gsf_v16");
  //leg_endcap->SetBorderSize(0);
  r_20_endcap->Draw();
  leg_endcap->Draw();
  r_30_endcap->Draw("same");
  r_40_endcap->Draw("same");
  r_50_endcap->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la_endcap->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la_endcap->DrawLatexNDC(0.6,0.65,"Endcap: 1.556 < |#eta| < 2.1");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("plots/debug_TriggerEfficiency_ele_%s_endcap.png", output_base.c_str()));


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
  rl1_7->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
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

  la2->DrawLatexNDC(0.22,0.92,"2025 Pythia8+Hydjet ZtoEE");
  la2->DrawLatexNDC(0.72,0.92,"L1 emulation");
  la2->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la2->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c2->SaveAs(Form("plots/debug_TriggerEfficiency_ele_%s_L1.png", output_base.c_str()));


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
  c3->SaveAs(Form("plots/debug_hiBin_ele_%s.png", output_base.c_str()));


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
  c4->SaveAs(Form("plots/debug_MissedEtaPhi_ele_%s.png", output_base.c_str()));


  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/debug_output_ele_%s.root", output_base.c_str()),"recreate");

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


