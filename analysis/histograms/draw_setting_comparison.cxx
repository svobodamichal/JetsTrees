// draw_setting_comparison.cxx

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"

#include <string>
#include <iostream>

using std::string;

static const int kNRadii = 3;
static const char* kRadii[kNRadii] = {"R0.2", "R0.3", "R0.4"};

static const int kNCentralities = 3;
static const char* kCentralities[kNCentralities] = {
  "CENT_0_10", "MID_20_40", "PERI_60_80"
};

static const int kNPtLeads = 4;
static const int kPtLeads[kNPtLeads] = {0, 5, 7, 9};

static const int kNSettings = 4;
static const char* kSettingTags[kNSettings] = {
  "MCRC1p5",
  "MCRC2p0",
  "MCVeto",
  "NoVeto"
};

static const char* kSettingDirs[kNSettings] = {
  "out_embedding_BAYES_Inclusive_MCRC1p5",
  "out_embedding_BAYES_Inclusive_MCRC2p0",
  "out_embedding_BAYES_Inclusive_MCVeto",
  "out_embedding_BAYES_Inclusive_NoVeto"
};

static TString NiceCentLabel(const string& centToken)
{
  int a = -1, b = -1;
  if (sscanf(centToken.c_str(), "CENT_%d_%d", &a, &b) == 2 ||
      sscanf(centToken.c_str(), "MID_%d_%d",  &a, &b) == 2 ||
      sscanf(centToken.c_str(), "PERI_%d_%d", &a, &b) == 2) {
    return Form("%d#font[52]{#minus}%d %%", a, b);
  }
  return TString(centToken.c_str());
}

static void StyleHist(TH1* h, int i)
{
  static const int cols[] = {kBlack, kRed+1, kAzure+2, kGreen+2};
  static const int mars[] = {20, 21, 22, 33};

  h->SetLineColor(cols[i]);
  h->SetMarkerColor(cols[i]);
  h->SetMarkerStyle(mars[i]);
  h->SetMarkerSize(1.1);
  h->SetLineWidth(2);
}

static TH1D* GetIter4(TFile* f,
                      const string& R,
                      const string& C,
                      int ptlead,
                      const char* setting)
{
  const string path = R + "_" + C + Form("_ptlead%d/Unfolded_iter4", ptlead);

  TH1D* h = dynamic_cast<TH1D*>(f->Get(path.c_str()));
  if (!h) {
    std::cout << "[missing] " << setting << " : " << path << std::endl;
    return 0;
  }

  TH1D* hc = dynamic_cast<TH1D*>(
    h->Clone(Form("h_%s_%s_%s_ptlead%d_iter4",
                  setting, R.c_str(), C.c_str(), ptlead))
  );

  hc->SetDirectory(0);
  return hc;
}

static void DrawOne(TFile* files[kNSettings],
                    const string& outDir,
                    const string& R,
                    const string& C,
                    int ptlead)
{
  TH1D* hs[kNSettings];

  for (int is = 0; is < kNSettings; ++is) {
    hs[is] = 0;
    if (!files[is]) continue;

    hs[is] = GetIter4(files[is], R, C, ptlead, kSettingTags[is]);
    if (!hs[is]) continue;

    StyleHist(hs[is], is);
    hs[is]->SetTitle("");
  }

  TH1D* hRef = hs[0];
  if (!hRef) {
    std::cout << "[skip] Missing reference for "
              << R << " " << C << " ptlead" << ptlead << std::endl;
    return;
  }

    // --------------------------------------------------
    // Optional shape normalization
    // false = absolute comparison, keeps yield differences
    // true  = normalize all spectra to the MCRC1p5 integral,
    //         so the ratio shows only pT-shape differences
    // --------------------------------------------------
    static const bool kNormalizeToReferenceIntegral = true;

    if (kNormalizeToReferenceIntegral) {
    const double refInt = hRef->Integral();

    for (int is = 1; is < kNSettings; ++is) {
        if (!hs[is]) continue;

        const double thisInt = hs[is]->Integral();
        if (thisInt > 0.0 && refInt > 0.0) {
        hs[is]->Scale(refInt / thisInt);
        }
    }
    }

  TCanvas* c = new TCanvas(
    Form("c_setting_%s_%s_ptlead%d", R.c_str(), C.c_str(), ptlead),
    "", 850, 950
  );

  TPad* pTop = new TPad("pTop", "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad("pBot", "", 0.0, 0.00, 1.0, 0.30);

  pTop->SetLeftMargin(0.13);
  pTop->SetRightMargin(0.04);
  pTop->SetTopMargin(0.05);
  pTop->SetBottomMargin(0.02);

  pBot->SetLeftMargin(0.13);
  pBot->SetRightMargin(0.04);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.30);

  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  gPad->SetLogy();

  double ymax = 0.0;
  for (int is = 0; is < kNSettings; ++is) {
    if (!hs[is]) continue;
    if (hs[is]->GetMaximum() > ymax) ymax = hs[is]->GetMaximum();
  }

  TLegend* leg = new TLegend(0.55, 0.62, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  bool first = true;
  for (int is = 0; is < kNSettings; ++is) {
    TH1D* h = hs[is];
    if (!h) continue;

    h->GetXaxis()->SetLabelSize(0);
    h->GetXaxis()->SetTitleSize(0);
    h->GetYaxis()->SetTitle("Unfolded yield, iter. 4");
    h->GetYaxis()->SetTitleOffset(1.25);

    h->SetMaximum(ymax * 5.0);
    h->SetMinimum(1e-12);

    if (first) {
      h->Draw("E1");
      first = false;
    } else {
      h->Draw("E1 SAME");
    }

    leg->AddEntry(h, kSettingTags[is], "lp");
  }

  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.17, 0.23, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
  lat.DrawLatex(0.17, 0.17, Form("%s, %s", R.c_str(), NiceCentLabel(C).Data()));
  lat.DrawLatex(0.17, 0.11,
                Form("#it{p}_{T}^{lead} #geq %d GeV/#it{c}", ptlead));

  pBot->cd();

  TH1D* firstRatio = 0;

  for (int is = 0; is < kNSettings; ++is) {
    if (!hs[is]) continue;

    TH1D* r = dynamic_cast<TH1D*>(
      hs[is]->Clone(Form("ratio_%s_%s_%s_ptlead%d",
                         kSettingTags[is], R.c_str(), C.c_str(), ptlead))
    );

    r->SetDirectory(0);
    r->Divide(hRef);
    StyleHist(r, is);

    r->SetTitle("");
    r->GetYaxis()->SetTitle("ratio to MCRC1p5");
    r->GetYaxis()->SetRangeUser(0.5, 1.5);
    r->GetYaxis()->SetTitleSize(0.085);
    r->GetYaxis()->SetLabelSize(0.080);
    r->GetYaxis()->SetTitleOffset(0.75);
    r->GetYaxis()->SetNdivisions(505);

    r->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    r->GetXaxis()->SetTitleSize(0.110);
    r->GetXaxis()->SetLabelSize(0.090);
    r->GetXaxis()->SetTitleOffset(1.05);

    if (!firstRatio) {
      firstRatio = r;
      r->Draw("E1");
    } else {
      r->Draw("E1 SAME");
    }
  }

  if (firstRatio) {
    const double xmin = firstRatio->GetXaxis()->GetXmin();
    const double xmax = firstRatio->GetXaxis()->GetXmax();

    TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(1);
    l1->Draw();

    TLine* l2 = new TLine(xmin, 1.1, xmax, 1.1);
    l2->SetLineColor(kGray+1);
    l2->SetLineStyle(2);
    l2->Draw();

    TLine* l3 = new TLine(xmin, 0.9, xmax, 0.9);
    l3->SetLineColor(kGray+1);
    l3->SetLineStyle(2);
    l3->Draw();
  }

  gSystem->mkdir(outDir.c_str(), true);

  const string outBase =
    outDir + "/setting_compare_iter4_" + R + "_" + C + Form("_ptlead%d", ptlead);

  c->SaveAs((outBase + ".pdf").c_str());
  c->SaveAs((outBase + ".png").c_str());

  // Intentionally no cleanup here.
  // ROOT 5 can segfault when canvases/pads own drawn objects during shutdown.
}

void draw_setting_comparison(const char* unfoldingDir,
                             const char* outDir = "setting_comparison_plots")
{
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* files[kNSettings];

  for (int is = 0; is < kNSettings; ++is) {
    files[is] = 0;

    const string path =
      string(unfoldingDir) + "/" + kSettingDirs[is] + "/responses_embedding.root";

    std::cout << "Opening " << kSettingTags[is] << " : " << path << std::endl;

    files[is] = TFile::Open(path.c_str(), "READ");
    if (!files[is] || files[is]->IsZombie()) {
      std::cout << "[error] Cannot open " << path << std::endl;
      files[is] = 0;
    }
  }

  if (!files[0]) {
    std::cout << "[fatal] Reference file MCRC1p5 missing. Stop." << std::endl;
    return;
  }

  for (int iR = 0; iR < kNRadii; ++iR) {
    for (int iC = 0; iC < kNCentralities; ++iC) {
      for (int ip = 0; ip < kNPtLeads; ++ip) {
        DrawOne(files, outDir, kRadii[iR], kCentralities[iC], kPtLeads[ip]);
      }
    }
  }

  std::cout << "Done drawing setting comparison plots." << std::endl;
}