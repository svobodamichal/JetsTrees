// draw_ptlead_comparison.cxx

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
#include <vector>
#include <iostream>

using std::string;
using std::vector;

static const int kNRadii = 3;
static const char* kRadii[kNRadii] = {
  "R0.2", "R0.3", "R0.4"
};

static const int kNCentralities = 3;
static const char* kCentralities[kNCentralities] = {
  "CENT_0_10", "MID_20_40", "PERI_60_80"
};

static const int kNPtLeads = 4;
static const int kPtLeads[kNPtLeads] = {
  0, 5, 7, 9
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

static TH1D* GetHist(TFile* f, const string& R, const string& C, int ptlead, const char* obj)
{
  const string path = R + "_" + C + Form("_ptlead%d/%s", ptlead, obj);
  TH1D* h = dynamic_cast<TH1D*>(f->Get(path.c_str()));

  if (!h) {
    std::cout << "[missing] " << path << std::endl;
    return 0;
  }

  TH1D* hc = dynamic_cast<TH1D*>(h->Clone(Form("%s_%s_%s_ptlead%d_clone", obj, R.c_str(), C.c_str(), ptlead)));
  hc->SetDirectory(0);
  return hc;
}

static void DrawOne(TFile* f,
                    const string& outDir,
                    const string& R,
                    const string& C,
                    const char* objName,
                    const char* label)
{
  vector<TH1D*> hs;
  vector<TH1D*> ratios;

  for (int i = 0; i < kNPtLeads; ++i) {
    TH1D* h = GetHist(f, R, C, kPtLeads[i], objName);
    if (!h) continue;

    StyleHist(h, i);
    h->SetTitle("");
    hs.push_back(h);
  }

  if (hs.empty()) return;

  TH1D* href = GetHist(f, R, C, 0, objName);
  if (!href) return;

  TCanvas* c = new TCanvas(Form("c_%s_%s_%s", objName, R.c_str(), C.c_str()), "", 850, 950);

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
  for (size_t ih = 0; ih < hs.size(); ++ih) {
    TH1D* h = hs[ih];
    if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
  }
  bool first = true;
  TLegend* leg = new TLegend(0.55, 0.62, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  for (size_t i = 0; i < hs.size(); ++i) {
    TH1D* h = hs[i];

    h->GetXaxis()->SetLabelSize(0);
    h->GetXaxis()->SetTitleSize(0);
    h->GetYaxis()->SetTitle(label);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->SetMaximum(ymax * 5.0);
    h->SetMinimum(1e-12);

    if (first) {
      h->Draw("E1");
      first = false;
    } else {
      h->Draw("E1 SAME");
    }

    leg->AddEntry(h, Form("#it{p}_{T}^{lead} #geq %d GeV/#it{c}", kPtLeads[i]), "lp");
  }

  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.17, 0.23, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
  lat.DrawLatex(0.17, 0.17, Form("%s, %s", R.c_str(), NiceCentLabel(C).Data()));
  lat.DrawLatex(0.17, 0.11, objName);

  pBot->cd();

  TH1D* firstRatio = 0;

  for (size_t i = 0; i < hs.size(); ++i) {
    TH1D* r = dynamic_cast<TH1D*>(hs[i]->Clone(Form("ratio_%d", (int)i)));
    r->SetDirectory(0);
    r->Divide(href);
    StyleHist(r, i);

    r->SetTitle("");
    r->GetYaxis()->SetTitle("ratio to #it{p}_{T}^{lead} #geq 0");
    r->GetYaxis()->SetRangeUser(0.0, 1.5);
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

    ratios.push_back(r);
  }

  if (firstRatio) {
    const double xmin = firstRatio->GetXaxis()->GetXmin();
    const double xmax = firstRatio->GetXaxis()->GetXmax();

    TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(1);
    l1->Draw();
  }

  gSystem->mkdir(outDir.c_str(), true);

  const string outBase = outDir + "/" + string(objName) + "_ptlead_compare_" + R + "_" + C;
  c->SaveAs((outBase + ".pdf").c_str());
  c->SaveAs((outBase + ".png").c_str());

//   delete c;
//   delete leg;
//   delete href;

// for (size_t ih = 0; ih < hs.size(); ++ih) delete hs[ih];
// for (size_t ir = 0; ir < ratios.size(); ++ir) delete ratios[ir];
}

void draw_ptlead_comparison(const char* inputRoot,
                            const char* outDir = "ptlead_comparison_plots")
{
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* f = TFile::Open(inputRoot, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "[error] Cannot open input file: " << inputRoot << std::endl;
    return;
  }

for (int iR = 0; iR < kNRadii; ++iR) {
  for (int iC = 0; iC < kNCentralities; ++iC) {
    DrawOne(f, outDir, kRadii[iR], kCentralities[iC], "hMeasFull", "Reco yield");
    DrawOne(f, outDir, kRadii[iR], kCentralities[iC], "hTrueFull", "MC truth yield");
    DrawOne(f, outDir, kRadii[iR], kCentralities[iC], "Unfolded_iter4", "Unfolded yield");
  }
}

  f->Close();
  delete f;
}