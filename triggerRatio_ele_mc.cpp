#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>

#include "pinchun/include/CommandLine.h"

namespace {

using HistPtr = std::unique_ptr<TH1D>;

HistPtr ReadHistogram(TFile &file, const char *name)
{
  TH1D *hist = dynamic_cast<TH1D*>(file.Get(name));
  if (!hist) {
    std::cerr << "Missing histogram " << name << " in " << file.GetName() << std::endl;
    return nullptr;
  }

  HistPtr clone(dynamic_cast<TH1D*>(hist->Clone()));
  clone->SetDirectory(nullptr);
  return clone;
}

void ApplyCurveStyle(TH1D &hist, Color_t color, Style_t markerStyle)
{
  hist.SetLineColor(color);
  hist.SetMarkerColor(color);
  hist.SetLineWidth(1.8);
  hist.SetMarkerSize(1.6);
  hist.SetMarkerStyle(markerStyle);
  hist.SetStats(0);
  hist.SetTitle("");
  hist.GetXaxis()->SetTitleSize(0.05);
  hist.GetYaxis()->SetTitleSize(0.05);
  hist.GetYaxis()->SetRangeUser(0.65, 1.35);
  hist.GetXaxis()->SetTitle("Reco electron #font[52]{p}_{T} [GeV]");
  hist.GetYaxis()->SetTitle("Efficiency ratio");
}

void DrawRatioPlot(const std::array<TH1D*, 5> &ratios,
                   const std::string &outputPath,
                   const std::string &nametag,
                   const std::string &ratioLabel,
                   const std::string &regionLabel,
                   float minHiBin,
                   float maxHiBin)
{
  TCanvas canvas("c_ratio", "c_ratio", 700, 600);
  TPad pad("p_ratio", "p_ratio", 0, 0, 1, 1);
  pad.SetLeftMargin(0.13);
  pad.SetBottomMargin(0.14);
  pad.Draw();
  pad.cd();

  TLegend legend(0.55, 0.165, 0.88, 0.445);
  legend.AddEntry(ratios[0], "HLT_HIEle15Gsf");
  legend.AddEntry(ratios[1], "HLT_HIEle20Gsf");
  legend.AddEntry(ratios[2], "HLT_HIEle30Gsf");
  legend.AddEntry(ratios[3], "HLT_HIEle40Gsf");
  legend.AddEntry(ratios[4], "HLT_HIEle50Gsf");

  ratios[0]->Draw("EP");
  TLine unityLine(ratios[0]->GetXaxis()->GetXmin(), 1.0, ratios[0]->GetXaxis()->GetXmax(), 1.0);
  unityLine.SetLineStyle(2);
  unityLine.SetLineWidth(2);
  unityLine.Draw("SAME");
  ratios[1]->Draw("EP SAME");
  ratios[2]->Draw("EP SAME");
  ratios[3]->Draw("EP SAME");
  ratios[4]->Draw("EP SAME");
  legend.Draw();

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.03);
  latex.DrawLatexNDC(0.18, 0.92, nametag.c_str());
  latex.DrawLatexNDC(0.18, 0.87, ratioLabel.c_str());
  latex.DrawLatexNDC(0.18, 0.82, regionLabel.c_str());
  latex.DrawLatexNDC(0.18, 0.77, Form("%.0f%% - %.0f%%", minHiBin / 2.0f, maxHiBin / 2.0f));

  canvas.SaveAs(outputPath.c_str());
}

void PrintUsage(const char *self)
{
  std::cout
    << "Usage: " << self << " [options]\n"
    << "  --inputNumerator <file>    ROOT output from triggerAnalysis for numerator menu\n"
    << "  --inputDenominator <file>  ROOT output from triggerAnalysis for denominator menu\n"
    << "  --outputBase <name>        Output tag for plots and ROOT file\n"
    << "  --nametag <label>          Plot label\n"
    << "  --ratioLabel <label>       Ratio label shown on plot\n"
    << "  --minHiBin <float>         Minimum hiBin\n"
    << "  --maxHiBin <float>         Maximum hiBin\n";
}

}  // namespace

int main(int argc, char *argv[])
{
  CommandLine cl(argc, argv);

  if (cl.GetBool("help", false) || cl.GetBool("h", false)) {
    PrintUsage(argv[0]);
    return 0;
  }

  const std::string inputNumerator = cl.Get("inputNumerator", "");
  const std::string inputDenominator = cl.Get("inputDenominator", "");
  const std::string outputBase = cl.Get("outputBase", "MCZee_ratio");
  const std::string nametag = cl.Get("nametag", "2025 Pythia8+Hydjet Z->EE");
  const std::string ratioLabel = cl.Get("ratioLabel", "HOverE0p5 / official");
  const float minHiBin = cl.GetDouble("minHiBin", 0.0);
  const float maxHiBin = cl.GetDouble("maxHiBin", 200.0);

  if (inputNumerator.empty() || inputDenominator.empty()) {
    std::cerr << "Both --inputNumerator and --inputDenominator are required." << std::endl;
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("output", kTRUE);

  TFile numeratorFile(inputNumerator.c_str(), "READ");
  TFile denominatorFile(inputDenominator.c_str(), "READ");
  if (numeratorFile.IsZombie() || denominatorFile.IsZombie()) {
    std::cerr << "Could not open one or both input ROOT files." << std::endl;
    return 1;
  }

  const std::array<const char*, 5> barrelNames = {"r_15", "r_20", "r_30", "r_40", "r_50"};
  const std::array<const char*, 5> endcapNames = {"r_15_endcap", "r_20_endcap", "r_30_endcap", "r_40_endcap", "r_50_endcap"};
  const std::array<Color_t, 5> colors = {kOrange + 7, kRed - 4, kBlue - 4, kGreen + 2, kMagenta - 9};
  const std::array<Style_t, 5> markers = {29, 20, 21, 22, 23};

  std::array<HistPtr, 5> barrelRatios;
  std::array<HistPtr, 5> endcapRatios;

  for (size_t i = 0; i < barrelNames.size(); ++i) {
    auto numeratorBarrel = ReadHistogram(numeratorFile, barrelNames[i]);
    auto denominatorBarrel = ReadHistogram(denominatorFile, barrelNames[i]);
    auto numeratorEndcap = ReadHistogram(numeratorFile, endcapNames[i]);
    auto denominatorEndcap = ReadHistogram(denominatorFile, endcapNames[i]);
    if (!numeratorBarrel || !denominatorBarrel || !numeratorEndcap || !denominatorEndcap) {
      return 1;
    }

    barrelRatios[i].reset(dynamic_cast<TH1D*>(numeratorBarrel->Clone(Form("%s_ratio", barrelNames[i]))));
    endcapRatios[i].reset(dynamic_cast<TH1D*>(numeratorEndcap->Clone(Form("%s_ratio", endcapNames[i]))));
    barrelRatios[i]->SetDirectory(nullptr);
    endcapRatios[i]->SetDirectory(nullptr);
    barrelRatios[i]->Divide(denominatorBarrel.get());
    endcapRatios[i]->Divide(denominatorEndcap.get());

    ApplyCurveStyle(*barrelRatios[i], colors[i], markers[i]);
    ApplyCurveStyle(*endcapRatios[i], colors[i], markers[i]);
  }

  DrawRatioPlot({barrelRatios[0].get(), barrelRatios[1].get(), barrelRatios[2].get(), barrelRatios[3].get(), barrelRatios[4].get()},
                Form("plots/TriggerEfficiencyRatio_ele_%s_barrel.png", outputBase.c_str()),
                nametag,
                ratioLabel,
                "Barrel: |#eta| < 1.442",
                minHiBin,
                maxHiBin);
  DrawRatioPlot({endcapRatios[0].get(), endcapRatios[1].get(), endcapRatios[2].get(), endcapRatios[3].get(), endcapRatios[4].get()},
                Form("plots/TriggerEfficiencyRatio_ele_%s_endcap.png", outputBase.c_str()),
                nametag,
                ratioLabel,
                "Endcap: 1.556 < |#eta| < 2.1",
                minHiBin,
                maxHiBin);

  auto outputFile = TFile::Open(Form("output/output_ele_ratio_%s.root", outputBase.c_str()), "RECREATE");
  for (const auto &hist : barrelRatios) hist->Write();
  for (const auto &hist : endcapRatios) hist->Write();
  outputFile->Close();

  return 0;
}
