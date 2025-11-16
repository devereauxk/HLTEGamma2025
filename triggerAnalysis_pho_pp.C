using namespace std;

Float_t EtMin = 10.;
Float_t EtMax = 100.;
Int_t NEtBins = 30;

Float_t etaMin = -5.0;
Float_t etaMax = 5.0;
Float_t NEtaBins = 100;

Float_t phiMin = -TMath::Pi();
Float_t phiMax = TMath::Pi();
Float_t NPhiBins = 100;

TH1D *denom = new TH1D("denom","denom",NEtBins,EtMin,EtMax);
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

TH1D *num_33 = new TH1D("num_33","num_33",NEtBins,EtMin,EtMax);
TH1D *numEta_40 = new TH1D("numEta_40","numEta_40",NEtaBins,etaMin,etaMax);
TH1D *numPhi_40 = new TH1D("numPhi_40","numPhi_40",NPhiBins,phiMin,phiMax);
TH1D *r_33 = new TH1D("r_33","r_33",NEtBins,EtMin,EtMax);

TH1D *num_60 = new TH1D("num_60","num_60",NEtBins,EtMin,EtMax);
TH1D *numEta_60 = new TH1D("numEta_60","numEta_60",NEtaBins,etaMin,etaMax);
TH1D *numPhi_60 = new TH1D("numPhi_60","numPhi_60",NPhiBins,phiMin,phiMax);
TH1D *r_60 = new TH1D("r_60","r_60",NEtBins,EtMin,EtMax);

TH1D *num_70 = new TH1D("num_70","num_70",NEtBins,EtMin,EtMax);
TH1D *numEta_70 = new TH1D("numEta_70","numEta_70",NEtaBins,etaMin,etaMax);
TH1D *numPhi_70 = new TH1D("numPhi_70","numPhi_70",NPhiBins,phiMin,phiMax);
TH1D *r_70 = new TH1D("r_70","r_70",NEtBins,EtMin,EtMax);

TH1D *num_80 = new TH1D("num_80","num_80",NEtBins,EtMin,EtMax);
TH1D *numEta_80 = new TH1D("numEta_80","numEta_80",NEtaBins,etaMin,etaMax);
TH1D *numPhi_80 = new TH1D("numPhi_80","numPhi_80",NPhiBins,phiMin,phiMax);
TH1D *r_80 = new TH1D("r_80","r_80",NEtBins,EtMin,EtMax);

TH1D *num_100 = new TH1D("num_100","num_100",NEtBins,EtMin,EtMax);
TH1D *numEta_100 = new TH1D("numEta_100","numEta_100",NEtaBins,etaMin,etaMax);
TH1D *numPhi_100 = new TH1D("numPhi_100","numPhi_100",NPhiBins,phiMin,phiMax);
TH1D *r_100 = new TH1D("r_100","r_100",NEtBins,EtMin,EtMax);

TH1D *num_120 = new TH1D("num_120","num_120",NEtBins,EtMin,EtMax);
TH1D *numEta_120 = new TH1D("numEta_120","numEta_120",NEtaBins,etaMin,etaMax);
TH1D *numPhi_120 = new TH1D("numPhi_120","numPhi_120",NPhiBins,phiMin,phiMax);
TH1D *r_120 = new TH1D("r_120","r_120",NEtBins,EtMin,EtMax);




void triggerAnalysis_pho_prompt(string input = "/eos/cms/store/group/phys_heavyions/jdlang/Run3_2025LowPUpp_ExpressForests/LowPUpp_SpecialHLTPhysics0_398683/crab_LowPUpp_SpecialHLTPhysics0_398683/251030_182029/0000/", string output_base = "lowPUpp_SpecialHLTPhysics0", int nfiles = -1){

  std::cout << "running triggerAnalysis_pho_prompt()" << std::endl;
  std::cout << "inputFile   = " << input  << std::endl;
  std::cout << "outputFile  = " << output_base << "_output.root" << std::endl;	

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");
  
  std::cout<< "Adding input files...";
  if (nfiles == -1) {
    string this_forest_name = input + "/HiForest_2025LowPUpp_*.root";
    HltTree   ->Add(this_forest_name.c_str());
    EventTree ->Add(this_forest_name.c_str());
    HiTree    ->Add(this_forest_name.c_str());
  }
  else {
    for(int i = 1; i <= nfiles; i++) {
      string this_forest_name = input + "/HiForest_2025LowPUpp_" + to_string(i) + ".root";
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


  // ================= HLT Tree =================
  Int_t           HLT_Photon33;
  Int_t           HLT_Photon50;
  Int_t           HLT_Photon75;
  Int_t           HLT_Photon33_PrescaleNumerator;
  Int_t           HLT_Photon33_PrescaleDenominator;

  HltTree->SetBranchStatus("*",0);     // disable all branches

  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  HltTree->SetBranchStatus("HLT_Photon33_v16", 1);
  HltTree->SetBranchStatus("HLT_Photon50_v24", 1);
  HltTree->SetBranchStatus("HLT_Photon75_v24", 1);

  HltTree->SetBranchAddress("HLT_Photon33_v16", &HLT_Photon33);
  HltTree->SetBranchAddress("HLT_Photon50_v24", &HLT_Photon50);
  HltTree->SetBranchAddress("HLT_Photon75_v24", &HLT_Photon75);

  HltTree->SetBranchAddress("HLT_Photon33_v16_PrescaleNumerator", &HLT_Photon33_PrescaleNumerator);
  HltTree->SetBranchAddress("HLT_Photon33_v16_PrescaleDenominator", &HLT_Photon33_PrescaleDenominator);

  // ================= Event Tree =================
  vector<float>   *phoEt = nullptr;
  vector<float>   *phoEta = nullptr;
  vector<float>   *phoPhi = nullptr;
  Int_t           nPho;
  vector<float>   *pfpIso3subUE = nullptr;
  vector<float>   *pfcIso3subUE = nullptr;
  vector<float>   *pfnIso3subUE = nullptr;
  vector<float>   *phoSCEta = nullptr;
  vector<float>   *phoSCPhi = nullptr;
  vector<float>   *phoHoverE = nullptr;
  vector<float>   *phoSigmaIEtaIEta_2012 = nullptr;

  //EventTree->Print();

  EventTree->SetBranchStatus("*",0);

  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

  EventTree->SetBranchStatus("nPho",1);
  //EventTree->SetBranchStatus("nEle",1);
  //EventTree->SetBranchStatus("mcPID",1);
  //EventTree->SetBranchStatus("pho_genMatchedIndex",1);
  //EventTree->SetBranchStatus("mcIndex",1);
  //EventTree->SetBranchStatus("mcEta",1);
  //EventTree->SetBranchStatus("mcPhi",1);

  EventTree->SetBranchStatus("phoEt",1);
  EventTree->SetBranchStatus("phoEta",1);
  EventTree->SetBranchStatus("phoHoverE",1);
  EventTree->SetBranchStatus("phoPhi",1);
  EventTree->SetBranchStatus("phoSCEta",1);
  EventTree->SetBranchStatus("phoSCPhi",1);

  EventTree->SetBranchStatus("phoSigmaIEtaIEta_2012",1);

  //EventTree->SetBranchStatus("mcCalIsoDR04",1);
  EventTree->SetBranchStatus("pfpIso3subUE",1);
  EventTree->SetBranchStatus("pfcIso3subUE",1);
  EventTree->SetBranchStatus("pfnIso3subUE",1);

  //EventTree->SetBranchStatus("pho_swissCrx",1);
  //EventTree->SetBranchStatus("pho_seedTime",1);

  //EventTree->SetBranchStatus("pho_ecalClusterIsoR3",1);
  //EventTree->SetBranchStatus("pho_hcalRechitIsoR3",1);
  //EventTree->SetBranchStatus("pho_trackIsoR3PtCut20",1);

  EventTree->SetBranchAddress("phoEt",&phoEt);
  EventTree->SetBranchAddress("phoEta",&phoEta);
  EventTree->SetBranchAddress("phoPhi",&phoPhi);
  EventTree->SetBranchAddress("nPho",&nPho);
  EventTree->SetBranchAddress("pfpIso3subUE",&pfpIso3subUE);
  EventTree->SetBranchAddress("pfcIso3subUE",&pfcIso3subUE);
  EventTree->SetBranchAddress("pfnIso3subUE",&pfnIso3subUE);
  EventTree->SetBranchAddress("phoSCEta",&phoSCEta);
  EventTree->SetBranchAddress("phoSCPhi",&phoSCPhi);
  EventTree->SetBranchAddress("phoHoverE",&phoHoverE);
  EventTree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);

  //ggHiNtuplizer ggHiNtuple;
  //ggHiNtuple.setupTreeForReading(EventTree);

  std::cout << "done" << std::endl;
  Long64_t entriesHLT = HltTree->GetEntries();
  std::cout << "HLT entries = " << entriesHLT << std::endl;
  Long64_t entriesTmp = HiTree->GetEntries();
  std::cout << "reco entries = " << entriesTmp << std::endl;


 

  // ================= Event Loop  =================

  // loop through HLT and create a key for each event
  //  std::map<unsigned long long, int> runLumiEvtToEntryMap; 

  //for(Long64_t i_entry = 0; i_entry < entriesHLT; i_entry++){

  // treeTrig->GetEntry(i_entry);
  // unsigned long long key = keyFromRunLumiEvent(hlt_run, hlt_lumi, hlt_event);
  // runLumiEvtToEntryMap[key] = i_entry;

  // }


  // loop through reco objects
  for (Long64_t i_event = 0; i_event < HltTree->GetEntries(); ++i_event){

    if(i_event%(HltTree->GetEntries()/500)==0) std::cout << "Processing entry " << i_event << " / " << entriesTmp << "\r" << std::flush;

    HiTree->GetEntry(i_event);
    HltTree->GetEntry(i_event);
    EventTree->GetEntry(i_event);
  	
    // event cuts
    if(fabs(vz)>15.0) continue;
    //if(hiBin>180) continue;	
       
    auto i_event_status = EventTree->GetEntry(i_event);
        
    //if(i_event_status < 0) {
    //    std::cout << Form("bad entry %lld in jet tree",i_event) << std::endl;
    //    continue;
    //}
	
    //unsigned long long key = keyFromRunLumiEvent(run,lumi,evt);

    //long long i_entry = -1;
    //Long64_t i_entry = i_event;

    //if(runLumiEvtToEntryMap.count(key) == 0) continue; // skip reco event if there is no HLT event match
    //else i_entry = runLumiEvtToEntryMap.at(key);

    //std::cout << "i_entry = " << i_entry << std::endl;


    // ============= find leading photon ==============

    float maxEt_denom = 0;
    float maxEta_denom = 0;
    float maxPhi_denom = 0;
    int i_leading = -1;
    
    for(Int_t i_track = 0; i_track < nPho; i_track++){
      if(phoEt->at(i_track) > maxEt_denom) { // find the leading phoEt in the event
        maxEt_denom = phoEt->at(i_track);
        i_leading = i_track;
      }
    }
    if (i_leading == -1) continue;
    maxEta_denom = phoEta->at(i_leading);
    maxPhi_denom = phoPhi->at(i_leading);

    // =========== apply track rejection ============
    float sumIso = pfpIso3subUE->at(i_leading) + pfcIso3subUE->at(i_leading) + pfnIso3subUE->at(i_leading);

    // barrel cut
    
    if (abs(phoSCEta->at(i_leading)) < 1.44) {
      if (sumIso > 11.697505)                              continue;
      if (phoHoverE->at(i_leading) > 0.247665)             continue;
      if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.012186) continue;
    }
    else continue; 

    // endcap cut
    else if (abs(phoSCEta->at(i_leading)) > 1.57 && abs(phoSCEta->at(i_leading)) < 2.1) {
      //if (sumIso > 20.911811)                              continue;
      //if (phoHoverE->at(i_leading) > 0.398866)             continue;
      //if (phoSigmaIEtaIEta_2012->at(i_leading) > 0.044998) continue;
    }
    else continue;
    

    //if (abs(phoSCEta->at(i_leading)) > 1.44 && !(abs(phoSCEta->at(i_leading)) > 1.57 && abs(phoSCEta->at(i_leading)) < 2.1)) continue;

    //if (abs(phoSCEta->at(i_leading)) > 1.44) continue;


    // ============== fill histograms ===============

    float prescale_weight = HLT_Photon33_PrescaleDenominator > 0 ? float(HLT_Photon33_PrescaleNumerator)/float(HLT_Photon33_PrescaleDenominator) : 1.0;

    // fill denominator histograms
    denom->Fill(maxEt_denom);

    if(HLT_Photon33) num_33->Fill(maxEt_denom, prescale_weight);
    /*
    if(HLT_Photon33) {// only fill the numerator if the trigger is on.

      float maxEt_num = 0;
      float maxEta_num = 0;
      float maxPhi_num = 0;
            
      // now fill the numerator
      for(Int_t i_track = 0; i_track < nPho; i_track++) {

        if(phoEt->at(i_track) > maxEt_num) { // find the leading phoEt in the event
          maxEt_num = phoEt->at(i_track);
          maxEta_num = phoEta->at(i_track);
          maxPhi_num = phoPhi->at(i_track);
        }

      }

      if(maxEt_num > 0) num_33->Fill(maxEt_num, weight);

    }
      */

  } // end of event loop

  r_33->Divide(num_33,denom,1,1,"B");

  r_33->SetLineColor(kRed-4);
  r_60->SetLineColor(kBlue-4);
  r_80->SetLineColor(kGreen+2);
  r_100->SetLineColor(kMagenta-9);
  r_120->SetLineColor(kPink+6);

  r_33->SetMarkerColor(kRed-4);
  r_60->SetMarkerColor(kBlue-4);
  r_80->SetMarkerColor(kGreen+2);
  r_100->SetMarkerColor(kMagenta-9);
  r_120->SetMarkerColor(kPink+6);

  double line_width = 1.8;
  r_33->SetLineWidth(line_width);
  r_60->SetLineWidth(line_width);
  r_80->SetLineWidth(line_width);
  r_100->SetLineWidth(line_width);
  r_120->SetLineWidth(line_width);

  double marker_size = 1.6;
  r_33->SetMarkerSize(marker_size);
  r_60->SetMarkerSize(marker_size);
  r_80->SetMarkerSize(marker_size);
  r_100->SetMarkerSize(marker_size);
  r_120->SetMarkerSize(marker_size);
    
  r_33->SetMarkerStyle(20);
  r_60->SetMarkerStyle(21);
  r_80->SetMarkerStyle(22);
  r_100->SetMarkerStyle(23);
  r_120->SetMarkerStyle(34);

  TCanvas *c1 = new TCanvas("c1","c1",700,600);
  c1->cd();
  TPad *p1 = new TPad("p1","p1",0,0,1,1);
  p1->SetLeftMargin(0.13);
  p1->SetBottomMargin(0.14);
  p1->Draw();
  p1->cd();
  r_33->SetTitle("");
  r_33->SetStats(0);
  r_33->GetXaxis()->SetTitleSize(0.05);
  r_33->GetYaxis()->SetTitleSize(0.05);
  r_33->GetXaxis()->SetTitle("photon #font[52]{E}_{T} [GeV]");
  r_33->GetYaxis()->SetTitle("Trigger efficiency");
  TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
  leg->AddEntry(r_33,"HLT_Photon33_v16");
  //leg->AddEntry(r_60,"HLT_PFJet60_v37");
  //leg->AddEntry(r_70,"");
  //leg->AddEntry(r_80,"HLT_PFJet80_L1SingleJet60_v1");
  //leg->AddEntry(r_100,"HLT_PFJet110_v16");
  //leg->AddEntry(r_120,"HLT_PFJet140_v35");
  //leg->SetBorderSize(0);
  r_33->Draw();
  leg->Draw();
  //r_33->Draw("same");
  //r_60->Draw("same");
  //r_80->Draw("same");
  //r_100->Draw("same");
  //r_120->Draw("same");

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03); 

  la->DrawLatexNDC(0.22,0.92,"2025 PPRef #sqrt{#it{s}} = 5.36 TeV");
  la->DrawLatexNDC(0.72,0.92,"Run 398683");
  //la->DrawLatexNDC(0.6,0.69,"2025 Run 3 MC");
  la->DrawLatexNDC(0.6,0.69,"|#eta^{#gamma}| < 1.44");

  c1->SaveAs(Form("TriggerEfficiency_pho_%s.png", output_base.c_str()));

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

  num_33->Write();
  num_60->Write();
  num_70->Write();
  num_80->Write();
  num_100->Write();
  num_120->Write();

  numEta_40->Write();
  numEta_60->Write();
  numEta_70->Write();
  numEta_80->Write();
  numEta_100->Write();
  numEta_120->Write();

  numPhi_40->Write();
  numPhi_60->Write();
  numPhi_70->Write();
  numPhi_80->Write();
  numPhi_100->Write();
  numPhi_120->Write();

  r_33->Write();
  r_60->Write();
  r_70->Write();
  r_80->Write();
  r_100->Write();
  r_120->Write();

  wf->Close();

}


