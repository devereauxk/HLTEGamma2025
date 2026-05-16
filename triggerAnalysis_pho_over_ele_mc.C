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

TH1D *denom        = new TH1D("denom",        "denom",        NPtBins,PtMin,PtMax);
TH1D *denom_endcap = new TH1D("denom_endcap", "denom_endcap", NPtBins,PtMin,PtMax);

TH1D *num_10        = new TH1D("num_10",        "num_10",        NPtBins,PtMin,PtMax);
TH1D *num_10_endcap = new TH1D("num_10_endcap", "num_10_endcap", NPtBins,PtMin,PtMax);
TH1D *r_pho10_ele        = new TH1D("r_pho10_ele",        "r_pho10_ele",        NPtBins,PtMin,PtMax);
TH1D *r_pho10_ele_endcap = new TH1D("r_pho10_ele_endcap", "r_pho10_ele_endcap", NPtBins,PtMin,PtMax);

TH1D *num_20        = new TH1D("num_20",        "num_20",        NPtBins,PtMin,PtMax);
TH1D *num_20_endcap = new TH1D("num_20_endcap", "num_20_endcap", NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele        = new TH1D("r_pho20_ele",        "r_pho20_ele",        NPtBins,PtMin,PtMax);
TH1D *r_pho20_ele_endcap = new TH1D("r_pho20_ele_endcap", "r_pho20_ele_endcap", NPtBins,PtMin,PtMax);

TH1D *num_30        = new TH1D("num_30",        "num_30",        NPtBins,PtMin,PtMax);
TH1D *num_30_endcap = new TH1D("num_30_endcap", "num_30_endcap", NPtBins,PtMin,PtMax);
TH1D *r_pho30_ele        = new TH1D("r_pho30_ele",        "r_pho30_ele",        NPtBins,PtMin,PtMax);
TH1D *r_pho30_ele_endcap = new TH1D("r_pho30_ele_endcap", "r_pho30_ele_endcap", NPtBins,PtMin,PtMax);

TH1D *num_40        = new TH1D("num_40",        "num_40",        NPtBins,PtMin,PtMax);
TH1D *num_40_endcap = new TH1D("num_40_endcap", "num_40_endcap", NPtBins,PtMin,PtMax);
TH1D *r_pho40_ele        = new TH1D("r_pho40_ele",        "r_pho40_ele",        NPtBins,PtMin,PtMax);
TH1D *r_pho40_ele_endcap = new TH1D("r_pho40_ele_endcap", "r_pho40_ele_endcap", NPtBins,PtMin,PtMax);

//TH1D *h_hiBin = new TH1D("h_hiBin","h_hiBin",250,0,250);


// ==========================================================
// main: photon-over-electron trigger efficiency in MC
// Denominator: leading quality electron (2023 PbPb cuts).
// Numerator:   photon triggers firing on that electron.
// No HLT object matching — photon triggers select on electrons
// without a dedicated dR requirement here.
// ==========================================================

void triggerAnalysis_pho_over_ele_mc(

    string inputForest = "/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/",
    string output_base = "MCphoOverEle0_100",
    string nametag = "2025 Pythia8+Hydjet Z->EE",

    int nfiles = -1,
    float minHiBin = 0.0,
    float maxHiBin = 200.0,
    int year = 2025,
    string plotDir = "plots"
  ){

  const string hlt_pho10 = (year == 2026) ? "HLT_HIGEDPhoton10_v17" : "HLT_HIGEDPhoton10_v16";
  const string hlt_pho20 = (year == 2026) ? "HLT_HIGEDPhoton20_v17" : "HLT_HIGEDPhoton20_v16";
  const string hlt_pho30 = (year == 2026) ? "HLT_HIGEDPhoton30_v17" : "HLT_HIGEDPhoton30_v16";
  const string hlt_pho40 = (year == 2026) ? "HLT_HIGEDPhoton40_v17" : "HLT_HIGEDPhoton40_v16";

  std::cout << "running triggerAnalysis_pho_over_ele_mc()" << std::endl;
  std::cout << "input forest directory = " << inputForest << std::endl;
  std::cout << "output tag             = " << output_base << std::endl;

  TChain *HltTree   = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree    = new TChain("hiEvtAnalyzer/HiTree");

  std::cout << "Adding input files...";

  if (nfiles == -1) {
    string p = inputForest + "/HiForestMiniAOD_*";
    HltTree  ->Add(p.c_str());
    EventTree->Add(p.c_str());
    HiTree   ->Add(p.c_str());
  }
  else {
    for(int i = 1; i <= nfiles; i++) {
      string p = inputForest + "/HiForestMiniAOD_" + to_string(i) + ".root";
      HltTree  ->Add(p.c_str());
      EventTree->Add(p.c_str());
      HiTree   ->Add(p.c_str());
    }
  }
  std::cout << "done" << std::endl;

  std::cout << "Setting Event, lumi, and run branchAddresses...";

  // ================= Hi Tree =================
  ULong64_t  evt;
  UInt_t     lumi;
  UInt_t     run;
  Int_t      hiBin;
  Float_t    hiHF;
  Float_t    pthat_weight;
  Float_t    vz;

  HiTree->SetBranchStatus("*",0);
  HiTree->SetBranchStatus("evt",1);
  HiTree->SetBranchStatus("lumi",1);
  HiTree->SetBranchStatus("run",1);
  HiTree->SetBranchStatus("hiBin",1);
  HiTree->SetBranchStatus("hiHF",1);
  HiTree->SetBranchStatus("weight",1);
  HiTree->SetBranchStatus("vz",1);

  HiTree->SetBranchAddress("evt",    &evt);
  HiTree->SetBranchAddress("lumi",   &lumi);
  HiTree->SetBranchAddress("run",    &run);
  HiTree->SetBranchAddress("hiBin",  &hiBin);
  HiTree->SetBranchAddress("hiHF",   &hiHF);
  HiTree->SetBranchAddress("weight", &pthat_weight);
  HiTree->SetBranchAddress("vz",     &vz);


  // ================= HLT Tree =================
  Int_t     L1_SingleEG7;
  Int_t     L1_SingleEG15;
  Int_t     L1_SingleEG21;
  Int_t     HLT_HIGEDPhoton10;
  Int_t     HLT_HIGEDPhoton20;
  Int_t     HLT_HIGEDPhoton30;
  Int_t     HLT_HIGEDPhoton40;
  ULong64_t forest_Event;
  Int_t     forest_LumiBlock, forest_Run;

  HltTree->SetBranchStatus("*",0);
  HltTree->SetBranchStatus("Event",1);
  HltTree->SetBranchStatus("LumiBlock",1);
  HltTree->SetBranchStatus("Run",1);
  HltTree->SetBranchStatus("L1_SingleEG7_BptxAND",1);
  HltTree->SetBranchStatus("L1_SingleEG15_BptxAND",1);
  HltTree->SetBranchStatus("L1_SingleEG21_BptxAND",1);
  HltTree->SetBranchStatus(hlt_pho10.c_str(),1);
  HltTree->SetBranchStatus(hlt_pho20.c_str(),1);
  HltTree->SetBranchStatus(hlt_pho30.c_str(),1);
  HltTree->SetBranchStatus(hlt_pho40.c_str(),1);

  HltTree->SetBranchAddress("Event",               &forest_Event);
  HltTree->SetBranchAddress("LumiBlock",           &forest_LumiBlock);
  HltTree->SetBranchAddress("Run",                 &forest_Run);
  HltTree->SetBranchAddress("L1_SingleEG7_BptxAND",  &L1_SingleEG7);
  HltTree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15);
  HltTree->SetBranchAddress("L1_SingleEG21_BptxAND", &L1_SingleEG21);
  HltTree->SetBranchAddress(hlt_pho10.c_str(), &HLT_HIGEDPhoton10);
  HltTree->SetBranchAddress(hlt_pho20.c_str(), &HLT_HIGEDPhoton20);
  HltTree->SetBranchAddress(hlt_pho30.c_str(), &HLT_HIGEDPhoton30);
  HltTree->SetBranchAddress(hlt_pho40.c_str(), &HLT_HIGEDPhoton40);


  // ================= Event Tree =================
  vector<float> *elePt               = nullptr;
  vector<float> *eleEta              = nullptr;
  vector<float> *elePhi              = nullptr;
  Int_t          nEle;
  vector<float> *eleHoverE           = nullptr;
  vector<float> *eleSigmaIEtaIEta_2012 = nullptr;
  vector<float> *eledEtaSeedAtVtx    = nullptr;
  vector<float> *eledPhiAtVtx        = nullptr;
  vector<float> *eleEoverPInv        = nullptr;
  vector<int>   *eleMissHits         = nullptr;
  vector<float> *eleD0               = nullptr;
  vector<float> *eleSCEta            = nullptr;
  vector<float> *eleSCRawEn          = nullptr;
  vector<float> *eleDz               = nullptr;

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);
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
  EventTree->SetBranchStatus("eleD0",1);
  EventTree->SetBranchStatus("eleSCEta",1);
  EventTree->SetBranchStatus("eleSCRawEn",1);
  EventTree->SetBranchStatus("eleDz",1);

  EventTree->SetBranchAddress("nEle",                   &nEle);
  EventTree->SetBranchAddress("elePt",                  &elePt);
  EventTree->SetBranchAddress("eleEta",                 &eleEta);
  EventTree->SetBranchAddress("eleHoverE",              &eleHoverE);
  EventTree->SetBranchAddress("elePhi",                 &elePhi);
  EventTree->SetBranchAddress("eleSigmaIEtaIEta_2012",  &eleSigmaIEtaIEta_2012);
  EventTree->SetBranchAddress("eledEtaSeedAtVtx",       &eledEtaSeedAtVtx);
  EventTree->SetBranchAddress("eledPhiAtVtx",           &eledPhiAtVtx);
  EventTree->SetBranchAddress("eleEoverPInv",           &eleEoverPInv);
  EventTree->SetBranchAddress("eleMissHits",            &eleMissHits);
  EventTree->SetBranchAddress("eleD0",                  &eleD0);
  EventTree->SetBranchAddress("eleSCEta",               &eleSCEta);
  EventTree->SetBranchAddress("eleSCRawEn",             &eleSCRawEn);
  EventTree->SetBranchAddress("eleDz",                  &eleDz);

  std::cout << "done" << std::endl;
  Long64_t entriesTmp = HltTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


  // ================= Event Loop =================

  int neles       = 0;
  int npass_probe = 0;
  int npass_numer = 0;

  for (ULong64_t i_event = 0; i_event < (ULong64_t)entriesTmp; ++i_event){

    if(i_event%(1+entriesTmp/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    HiTree   ->GetEntry(i_event);
    HltTree  ->GetEntry(i_event);
    EventTree->GetEntry(i_event);

    if(fabs(vz) > 15.0) continue;
    if(hiBin < minHiBin || hiBin >= maxHiBin) continue;

    //h_hiBin->Fill(hiBin);
    neles += nEle;

    float weight = pthat_weight;

    // ====== single pass: find leading quality electron per region ======

    int   i_lead_barrel = -1;  float lead_barrel_pt = -1.f;
    int   i_lead_endcap = -1;  float lead_endcap_pt = -1.f;

    for(Int_t i_ele = 0; i_ele < nEle; i_ele++){

      bool isBarrel = false;
      bool isEndcap = false;

      // WP95 cut-based ID (max threshold across centrality bins)
      if (abs(eleEta->at(i_ele)) < 1.442) {
        if (eleMissHits->at(i_ele) > 2)                  continue;
        if (abs(eleD0->at(i_ele)) > 0.05)                continue;
        if (abs(eleDz->at(i_ele)) > 0.10)                continue;
        if (eleSigmaIEtaIEta_2012->at(i_ele) > 0.0131)   continue;
        if (abs(eledEtaSeedAtVtx->at(i_ele)) > 0.00457)  continue;
        if (abs(eledPhiAtVtx->at(i_ele)) > 0.0963)       continue;
        if (eleEoverPInv->at(i_ele) > 0.421)             continue;
        if (eleHoverE->at(i_ele) > 0.156)                continue;
        isBarrel = true;
      }
      else if (abs(eleEta->at(i_ele)) > 1.556 && abs(eleEta->at(i_ele)) < 2.1) {
        if (eleMissHits->at(i_ele) > 3)                  continue;
        if (abs(eleD0->at(i_ele)) > 0.10)                continue;
        if (abs(eleDz->at(i_ele)) > 0.20)                continue;
        if (eleSigmaIEtaIEta_2012->at(i_ele) > 0.0382)   continue;
        if (abs(eledEtaSeedAtVtx->at(i_ele)) > 0.00881)  continue;
        if (abs(eledPhiAtVtx->at(i_ele)) > 0.264)        continue;
        if (eleEoverPInv->at(i_ele) > 0.146)             continue;
        if (eleHoverE->at(i_ele) > 0.178)                continue;
        isEndcap = true;
      }
      else continue;

      if (isBarrel && elePt->at(i_ele) > lead_barrel_pt) {
        lead_barrel_pt = elePt->at(i_ele); i_lead_barrel = i_ele;
      }
      if (isEndcap && elePt->at(i_ele) > lead_endcap_pt) {
        lead_endcap_pt = elePt->at(i_ele); i_lead_endcap = i_ele;
      }

    } // end of per-electron loop

    if (i_lead_barrel == -1 && i_lead_endcap == -1) continue;
    npass_probe++;

    // fill barrel (at most once per event)
    if (i_lead_barrel != -1) {
      float recoPt = elePt->at(i_lead_barrel);
      denom->Fill(recoPt, weight);

      if(HLT_HIGEDPhoton10) { num_10->Fill(recoPt, weight); npass_numer++; }
      if(HLT_HIGEDPhoton20)   num_20->Fill(recoPt, weight);
      if(HLT_HIGEDPhoton30)   num_30->Fill(recoPt, weight);
      if(HLT_HIGEDPhoton40)   num_40->Fill(recoPt, weight);
    }

    // fill endcap (at most once per event)
    if (i_lead_endcap != -1) {
      float recoPt = elePt->at(i_lead_endcap);
      denom_endcap->Fill(recoPt, weight);

      if(HLT_HIGEDPhoton10) num_10_endcap->Fill(recoPt, weight);
      if(HLT_HIGEDPhoton20) num_20_endcap->Fill(recoPt, weight);
      if(HLT_HIGEDPhoton30) num_30_endcap->Fill(recoPt, weight);
      if(HLT_HIGEDPhoton40) num_40_endcap->Fill(recoPt, weight);
    }

  } // end of event loop

  std::cout << std::endl;
  std::cout << "HLT_HIGEDPhoton10 results ..." << std::endl;
  std::cout << "neles       = " << neles       << std::endl;
  std::cout << "npass_probe = " << npass_probe << std::endl;
  std::cout << "npass_numer = " << npass_numer << std::endl;

  r_pho10_ele->Divide(num_10, denom, 1,1,"B");
  r_pho20_ele->Divide(num_20, denom, 1,1,"B");
  r_pho30_ele->Divide(num_30, denom, 1,1,"B");
  r_pho40_ele->Divide(num_40, denom, 1,1,"B");

  r_pho10_ele_endcap->Divide(num_10_endcap, denom_endcap, 1,1,"B");
  r_pho20_ele_endcap->Divide(num_20_endcap, denom_endcap, 1,1,"B");
  r_pho30_ele_endcap->Divide(num_30_endcap, denom_endcap, 1,1,"B");
  r_pho40_ele_endcap->Divide(num_40_endcap, denom_endcap, 1,1,"B");

  r_pho10_ele->SetLineColor(kRed-4);
  r_pho20_ele->SetLineColor(kBlue-4);
  r_pho30_ele->SetLineColor(kGreen+2);
  r_pho40_ele->SetLineColor(kMagenta-9);
  r_pho10_ele_endcap->SetLineColor(kRed-4);
  r_pho20_ele_endcap->SetLineColor(kBlue-4);
  r_pho30_ele_endcap->SetLineColor(kGreen+2);
  r_pho40_ele_endcap->SetLineColor(kMagenta-9);

  r_pho10_ele->SetMarkerColor(kRed-4);
  r_pho20_ele->SetMarkerColor(kBlue-4);
  r_pho30_ele->SetMarkerColor(kGreen+2);
  r_pho40_ele->SetMarkerColor(kMagenta-9);
  r_pho10_ele_endcap->SetMarkerColor(kRed-4);
  r_pho20_ele_endcap->SetMarkerColor(kBlue-4);
  r_pho30_ele_endcap->SetMarkerColor(kGreen+2);
  r_pho40_ele_endcap->SetMarkerColor(kMagenta-9);

  double line_width = 1.8;
  r_pho10_ele->SetLineWidth(line_width);
  r_pho20_ele->SetLineWidth(line_width);
  r_pho30_ele->SetLineWidth(line_width);
  r_pho40_ele->SetLineWidth(line_width);
  r_pho10_ele_endcap->SetLineWidth(line_width);
  r_pho20_ele_endcap->SetLineWidth(line_width);
  r_pho30_ele_endcap->SetLineWidth(line_width);
  r_pho40_ele_endcap->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_pho10_ele->SetMarkerSize(marker_size);
  r_pho20_ele->SetMarkerSize(marker_size);
  r_pho30_ele->SetMarkerSize(marker_size);
  r_pho40_ele->SetMarkerSize(marker_size);
  r_pho10_ele_endcap->SetMarkerSize(marker_size);
  r_pho20_ele_endcap->SetMarkerSize(marker_size);
  r_pho30_ele_endcap->SetMarkerSize(marker_size);
  r_pho40_ele_endcap->SetMarkerSize(marker_size);

  r_pho10_ele->SetMarkerStyle(20);
  r_pho20_ele->SetMarkerStyle(21);
  r_pho30_ele->SetMarkerStyle(22);
  r_pho40_ele->SetMarkerStyle(23);
  r_pho10_ele_endcap->SetMarkerStyle(20);
  r_pho20_ele_endcap->SetMarkerStyle(21);
  r_pho30_ele_endcap->SetMarkerStyle(22);
  r_pho40_ele_endcap->SetMarkerStyle(23);

  // =============== Fig HLT BARREL ==================

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
  r_pho10_ele->GetYaxis()->SetRangeUser(0,1.05);
  r_pho10_ele->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_pho10_ele->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_pho10_ele,"HLT_HIGEDPhoton10");
  leg->AddEntry(r_pho20_ele,"HLT_HIGEDPhoton20");
  leg->AddEntry(r_pho30_ele,"HLT_HIGEDPhoton30");
  leg->AddEntry(r_pho40_ele,"HLT_HIGEDPhoton40");
  r_pho10_ele->Draw();
  leg->Draw();
  r_pho20_ele->Draw("same");
  r_pho30_ele->Draw("same");
  r_pho40_ele->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->DrawLatexNDC(0.22,0.92,Form("%s", nametag.c_str()));
  la->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la->DrawLatexNDC(0.6,0.65,"Barrel: |#eta| < 1.442");
  la->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1->SaveAs(Form("%s/TriggerEfficiency_phoOverEle_%s_barrel.png", plotDir.c_str(), output_base.c_str()));


  // =============== Fig HLT ENDCAP ==================

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
  r_pho10_ele_endcap->GetYaxis()->SetRangeUser(0,1.05);
  r_pho10_ele_endcap->GetXaxis()->SetTitle("electron #font[52]{p}_{T} [GeV]");
  r_pho10_ele_endcap->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg_endcap = new TLegend(0.55,0.3,0.88,0.5);
  leg_endcap->AddEntry(r_pho10_ele_endcap,"HLT_HIGEDPhoton10");
  leg_endcap->AddEntry(r_pho20_ele_endcap,"HLT_HIGEDPhoton20");
  leg_endcap->AddEntry(r_pho30_ele_endcap,"HLT_HIGEDPhoton30");
  leg_endcap->AddEntry(r_pho40_ele_endcap,"HLT_HIGEDPhoton40");
  r_pho10_ele_endcap->Draw();
  leg_endcap->Draw();
  r_pho20_ele_endcap->Draw("same");
  r_pho30_ele_endcap->Draw("same");
  r_pho40_ele_endcap->Draw("same");

  TLatex *la_endcap = new TLatex();
  la_endcap->SetTextFont(42);
  la_endcap->SetTextSize(0.03);
  la_endcap->DrawLatexNDC(0.22,0.92,Form("%s", nametag.c_str()));
  la_endcap->DrawLatexNDC(0.72,0.92,"no L1 emulation");
  la_endcap->DrawLatexNDC(0.6,0.65,"Endcap: 1.556 < |#eta| < 2.1");
  la_endcap->DrawLatexNDC(0.6,0.6,Form("%.0f%% - %.0f%%", minHiBin/2, maxHiBin/2));

  c1_endcap->SaveAs(Form("%s/TriggerEfficiency_phoOverEle_%s_endcap.png", plotDir.c_str(), output_base.c_str()));


  // ====== hiBin histogram — commented out ======
  //TCanvas *c3 = new TCanvas("c3","c3",700,600);
  //c3->cd();
  //TPad *p3 = new TPad("p3","p3",0,0,1,1);
  //p3->SetLeftMargin(0.13); p3->SetBottomMargin(0.14); p3->Draw(); p3->cd();
  //h_hiBin->SetTitle(""); h_hiBin->SetStats(0);
  //h_hiBin->GetXaxis()->SetTitle("hiBin"); h_hiBin->GetYaxis()->SetTitle("Counts");
  //h_hiBin->Draw();
  //c3->SaveAs(Form("%s/hiBin_phoOverEle_%s.png", plotDir.c_str(), output_base.c_str()));


  // =============== Save histograms ==================

  auto wf = TFile::Open(Form("output/output_phoOverEle_%s.root", output_base.c_str()),"recreate");

  denom->Write();
  denom_endcap->Write();

  num_10->Write(); num_20->Write(); num_30->Write(); num_40->Write();
  num_10_endcap->Write(); num_20_endcap->Write(); num_30_endcap->Write(); num_40_endcap->Write();

  r_pho10_ele->Write(); r_pho20_ele->Write(); r_pho30_ele->Write(); r_pho40_ele->Write();
  r_pho10_ele_endcap->Write(); r_pho20_ele_endcap->Write(); r_pho30_ele_endcap->Write(); r_pho40_ele_endcap->Write();

  //h_hiBin->Write();

  wf->Close();

}
