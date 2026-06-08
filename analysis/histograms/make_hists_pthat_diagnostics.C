#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <cmath>

// =========================================================
//  Jet-quality cuts; keep synchronized with analysis/maker
// =========================================================
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4
const double CUT_NEUTRAL_FRACTION = 0.95;
const double RECO_DUMMY_CUT = -500.0;

// Use the same pTlead thresholds as in the analysis code.
const int    N_LEAD = 4;
const double PTLEAD_THR[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

// pThat bin boundaries: [3,5), [5,7), ... [50,inf)
// In this tree we identify the pThat bin by xsecWeight, not by pThat itself.
static const int kNPthatBins = 11;
static const double kPtHatLow[kNPthatBins]  = {3, 5, 7, 9, 11, 15, 20, 25, 30, 40, 50};
static const double kPtHatHigh[kNPthatBins] = {5, 7, 9, 11, 15, 20, 25, 30, 40, 50, -1};

static const double kXsecWeights[kNPthatBins] = {
  1.616e+0, 1.355e-01, 2.288e-02, 5.524e-03, 2.203e-03,
  3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09
};

static const double kMinSignif = std::sqrt(10.0);

// Analysis binning: reco/corrected-reco and truth.
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static const int nbins_overlay = 60;
static const double xmin_overlay = 0.0;
static const double xmax_overlay = 60.0;

static double areaMinForR(double R)
{
  if (R < 0.25) return CUT_AREA_02;
  if (R < 0.35) return CUT_AREA_03;
  return CUT_AREA_04;
}

static int FindPtHatBin(double xsecW)
{
  const double relTol = 1e-6;
  for (int i = 0; i < kNPthatBins; ++i) {
    const double ref = kXsecWeights[i];
    const double diff = std::fabs(xsecW - ref);
    if (diff <= relTol * std::fabs(ref)) return i;
  }
  return -1;
}

static TString PtHatLabel(int ip)
{
  TString s;
  if (kPtHatHigh[ip] < 0) s.Form("%.0f < #hat{p}_{T}", kPtHatLow[ip]);
  else s.Form("%.0f < #hat{p}_{T} < %.0f", kPtHatLow[ip], kPtHatHigh[ip]);
  return s;
}

static TString PtHatName(int ip)
{
  TString s;
  if (kPtHatHigh[ip] < 0) s.Form("pthat%.0f_inf", kPtHatLow[ip]);
  else s.Form("pthat%.0f_%.0f", kPtHatLow[ip], kPtHatHigh[ip]);
  return s;
}

static void SetPthatStyle(TH1D* h, int ip)
{
  static const int colors[kNPthatBins] = {
    kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1,
    kCyan+2, kOrange+7, kViolet+7, kAzure+7, kPink+7, kSpring+5
  };
  static const int markers[kNPthatBins] = {
    20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 33
  };
  h->SetLineColor(colors[ip]);
  h->SetMarkerColor(colors[ip]);
  h->SetMarkerStyle(markers[ip]);
  h->SetMarkerSize(0.8);
  h->SetLineWidth(2);
}

static void DrawLatexBlock(double Rval, const TString& centLabel)
{
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.036);
  tex.DrawLatex(0.20, 0.24, Form("R = %.1f, %s", Rval, centLabel.Data()));
  tex.DrawLatex(0.20, 0.19, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
  tex.DrawLatex(0.20, 0.14, "embedding");
}

static void DrawOverlayPthat(TH1D* h[kNPthatBins],
                             const char* quantityLabel,
                             const char* xaxisTitle,
                             const TString& pdfName,
                             double Rval,
                             const TString& centLabel,
                             double xmin,
                             double xmax)
{
  TCanvas* c = new TCanvas("c_overlay_pthat", "c_overlay_pthat", 900, 750);
  c->SetLogy();
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);

  double ymax = 0.0;
  double ymin = 1e99;
  for (int ip = 0; ip < kNPthatBins; ++ip) {
    for (int ib = 1; ib <= h[ip]->GetNbinsX(); ++ib) {
      const double y = h[ip]->GetBinContent(ib);
      if (y > ymax) ymax = y;
      if (y > 0.0 && y < ymin) ymin = y;
    }
  }
  if (ymax <= 0.0) ymax = 1.0;
  if (ymin == 1e99) ymin = 1e-12;

  bool first = true;
  TLegend* leg = new TLegend(0.68, 0.46, 0.92, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.018);
  leg->SetNColumns(1);

  for (int ip = 0; ip < kNPthatBins; ++ip) {
    SetPthatStyle(h[ip], ip);
    h[ip]->SetTitle("");
    h[ip]->GetXaxis()->SetTitle(xaxisTitle);
    h[ip]->GetYaxis()->SetTitle("centrality-weighted counts");
    h[ip]->GetXaxis()->SetRangeUser(xmin, xmax);
    h[ip]->SetMinimum(ymin * 0.3);
    h[ip]->SetMaximum(ymax * 20.0);

    if (first) {
      h[ip]->Draw("E1");
      first = false;
    } else {
      h[ip]->Draw("E1 SAME");
    }
    leg->AddEntry(h[ip], PtHatLabel(ip), "lep");
  }

  leg->Draw();

  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextSize(0.036);
  title.DrawLatex(0.16, 0.92, quantityLabel);

 // DrawLatexBlock(Rval, centLabel);

  c->SaveAs(pdfName.Data());
  TString pngName = pdfName;
  pngName.ReplaceAll(".pdf", ".png");
  c->SaveAs(pngName.Data());

  delete leg;
  delete c;
}

static void DrawSignificanceMaskPlot(TH1D* hRaw,
                                     TH1D* hKeep,
                                     TH1D* hReject,
                                     TH1D* hEntries,
                                     int ip,
                                     const TString& pdfName,
                                     double Rval,
                                     const TString& centLabel,
                                     double ptleadThr)
{
double nAll  = 0.0;
double nBlue = 0.0;

for (int ix = 1; ix <= hRaw->GetNbinsX(); ++ix) {
  const double c = hRaw->GetBinContent(ix);
  const double e = hRaw->GetBinError(ix);
  const double signif = (e > 0.0 ? c / e : 0.0);
  const bool keep = (signif > kMinSignif);

  const double n = hEntries->GetBinContent(ix);

  nAll += n;
  if (keep) nBlue += n;
}


  TCanvas* cSig = new TCanvas("c_signif_mask", "c_signif_mask", 900, 750);
  cSig->SetLogy();
  cSig->SetLeftMargin(0.13);
  cSig->SetRightMargin(0.04);
  cSig->SetBottomMargin(0.12);

  hRaw->SetTitle("");
  hRaw->GetXaxis()->SetTitle("#it{p}_{T,reco}^{corr} (GeV/#it{c})");
  hRaw->GetYaxis()->SetTitle("centrality-weighted counts");
  hRaw->GetXaxis()->SetRangeUser(-40, 60);

  double ymax = hRaw->GetMaximum();
  if (ymax <= 0.0) ymax = 1.0;
  hRaw->SetMinimum(1e-12);
  hRaw->SetMaximum(ymax * 30.0);

  hRaw->SetLineColor(kGray+2);
  hRaw->SetMarkerColor(kGray+2);
  hRaw->SetMarkerStyle(24);
  hRaw->SetLineWidth(1);

  hKeep->SetLineColor(kBlue+1);
  hKeep->SetMarkerColor(kBlue+1);
  hKeep->SetMarkerStyle(20);
  hKeep->SetLineWidth(2);

  hReject->SetLineColor(kRed+1);
  hReject->SetMarkerColor(kRed+1);
  hReject->SetMarkerStyle(25);
  hReject->SetLineWidth(2);

  hRaw->Draw("E1");
  hKeep->Draw("E1 SAME");
  hReject->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.16, 0.70, 0.56, 0.88);
  leg->SetTextSize(0.026);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hRaw,
                Form("all bins: N_{jets} = %.0f", nAll),
                "lep");

  leg->AddEntry(hKeep,
                Form("accepted: N_{jets} = %.0f / %.0f", nBlue, nAll),
                "lep");

  leg->AddEntry(hReject,
                Form("rejected: N_{jets} = %.0f", nAll - nBlue),
                "lep");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.036);
  tex.DrawLatex(0.16, 0.92, Form("Significance mask, %s", PtHatLabel(ip).Data()));
  tex.DrawLatex(0.20, 0.25, Form("#it{p}_{T}^{lead} #geq %.0f GeV", ptleadThr));
  tex.DrawLatex(0.20, 0.20, Form("R = %.1f, %s", Rval, centLabel.Data()));

  cSig->SaveAs(pdfName.Data());
  TString pngName = pdfName;
  pngName.ReplaceAll(".pdf", ".png");
  cSig->SaveAs(pngName.Data());

  delete leg;
  delete cSig;
}

void make_hists_pthat_diagnostics(const char *infile  = "embedding_merged.root",
                                  const char *outfile = "hists_pthat_diagnostics.root")
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  const char* pdfTop = "pdf_pthat_diagnostics";
  gSystem->mkdir(pdfTop, kTRUE);

  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << std::endl;
    return;
  }


  TFile *fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create output file " << outfile << std::endl;
    fin->Close();
    return;
  }

  std::cout << "Input : " << infile  << std::endl;
  std::cout << "Output: " << outfile << std::endl;

  TIter nextR(fin->GetListOfKeys());
  TKey *keyR = 0;

  while ((keyR = (TKey*) nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) continue;

    const double dRmax = 0.6 * R;
    const double areaMin = areaMinForR(R);

    TDirectoryFile *dirR = (TDirectoryFile*) keyR->ReadObj();
    if (!dirR) continue;

    std::cout << "Radius directory: " << rname
              << " (R=" << R << ", dRmax=" << dRmax
              << ", areaMin=" << areaMin << ")\n";

    fout->cd();
    TDirectory *outR = fout->mkdir(rname.c_str());
    if (!outR) outR = fout->GetDirectory(rname.c_str());

    TString pdfRDir; pdfRDir.Form("%s/%s", pdfTop, rname.c_str());
    gSystem->mkdir(pdfRDir.Data(), kTRUE);

    TIter nextC(dirR->GetListOfKeys());
    TKey *keyC = 0;

    while ((keyC = (TKey*) nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      TDirectoryFile *dirC = (TDirectoryFile*) keyC->ReadObj();
      if (!dirC) continue;

      TTree *tree = (TTree*) dirC->Get("JetTree");
      if (!tree) continue;

      const Long64_t nentries = tree->GetEntries();
      std::cout << "  Centrality: " << cname << " (entries: " << nentries << ")\n";

      if (!tree->GetBranch("mc_pt") ||
          !tree->GetBranch("xsecWeight") ||
          !tree->GetBranch("deltaR") ||
          !tree->GetBranch("centralityWeight") ||
          !tree->GetBranch("reco_pt") ||
          !tree->GetBranch("reco_pt_corr") ||
          !tree->GetBranch("reco_area") ||
          !tree->GetBranch("reco_neutral_fraction") ||
          !tree->GetBranch("reco_trigger_match") ||
          !tree->GetBranch("reco_pt_lead")) {
        std::cerr << "    Missing branches in " << rname << "/" << cname << " -> skipping.\n";
        continue;
      }

      outR->cd();
      TDirectory *outC = outR->mkdir(cname.c_str());
      if (!outC) outC = outR->GetDirectory(cname.c_str());
      outC->cd();

      TDirectory *outPthat = outC->mkdir("pthat_diagnostics");
      if (!outPthat) outPthat = outC->GetDirectory("pthat_diagnostics");

      TString pdfCDir; pdfCDir.Form("%s/%s", pdfRDir.Data(), cname.c_str());
      gSystem->mkdir(pdfCDir.Data(), kTRUE);

      Float_t xsecWeight = 0.f;
      Float_t centralityWeight = 1.f;
      Float_t mc_pt = -999.f;
      Float_t deltaR = -1.f;
      Float_t reco_pt = -999.f;
      Float_t reco_pt_corr = -999.f;
      Float_t reco_area = -999.f;
      Float_t reco_neutral_fraction = 999.f;
      Float_t reco_pt_lead = -999.f;
      Bool_t  reco_trigger_match = kFALSE;

      tree->SetBranchAddress("xsecWeight",        &xsecWeight);
      tree->SetBranchAddress("centralityWeight",  &centralityWeight);
      tree->SetBranchAddress("mc_pt",             &mc_pt);
      tree->SetBranchAddress("deltaR",            &deltaR);
      tree->SetBranchAddress("reco_pt",           &reco_pt);
      tree->SetBranchAddress("reco_pt_corr",      &reco_pt_corr);
      tree->SetBranchAddress("reco_area",         &reco_area);
      tree->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);
      tree->SetBranchAddress("reco_pt_lead",      &reco_pt_lead);
      tree->SetBranchAddress("reco_trigger_match",&reco_trigger_match);

      // One set per pThat bin. These are intentionally filled with centralityWeight only.
      // The xsecWeight is used only to identify the bin. This lets you inspect per-bin shapes
      // before final xsec-weighted merging.
      TH1D* hMc_pthat[kNPthatBins];
      TH1D* hReco_pthat[kNPthatBins];

      // fine-binned overlay only
      TH1D* hRecoCorr_overlay[N_LEAD][kNPthatBins];

      // analysis-binned significance-mask histograms
      TH1D* hRecoCorr_pthat[N_LEAD][kNPthatBins];
      TH1D* hRecoCorr_keep[N_LEAD][kNPthatBins];
      TH1D* hRecoCorr_reject[N_LEAD][kNPthatBins];
      TH1D* hRecoCorr_entries[N_LEAD][kNPthatBins];

      for (int ip = 0; ip < kNPthatBins; ++ip) {
        hMc_pthat[ip] = new TH1D(Form("hMcPt_%s_%s_%s", PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                Form("MC p_{T}, %s;p_{T}^{MC} [GeV/#it{c}];centrality-weighted counts", PtHatLabel(ip).Data()),
                                nbins_overlay, xmin_overlay, xmax_overlay);

        hReco_pthat[ip] = new TH1D(Form("hRecoPt_%s_%s_%s", PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                  Form("Reco p_{T}, %s;p_{T}^{reco} [GeV/#it{c}];centrality-weighted counts", PtHatLabel(ip).Data()),
                                  nbins_overlay, xmin_overlay, xmax_overlay);

        for (int it = 0; it < N_LEAD; ++it) {
          hRecoCorr_overlay[it][ip] = new TH1D(Form("hRecoPtCorr_OVERLAY_ptlead%.0f_%s_%s_%s", PTLEAD_THR[it], PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                  Form("Reco p_{T}^{corr}, %s, p_{T}^{lead} >= %.0f;p_{T}^{reco,corr} [GeV/#it{c}];centrality-weighted counts", PtHatLabel(ip).Data(), PTLEAD_THR[it]),
                                  nbins_overlay, xmin_overlay, xmax_overlay);
          hRecoCorr_pthat[it][ip] = new TH1D(Form("hRecoPtCorr_ptlead%.0f_%s_%s_%s", PTLEAD_THR[it], PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                            Form("Reco p_{T}^{corr}, %s, p_{T}^{lead} >= %.0f;p_{T}^{reco,corr} [GeV/#it{c}];centrality-weighted counts",
                                                 PtHatLabel(ip).Data(), PTLEAD_THR[it]),
                                            nbins_meas, bin_meas_edges);

          hRecoCorr_keep[it][ip] = new TH1D(Form("hRecoPtCorr_ACCEPTED_ptlead%.0f_%s_%s_%s", PTLEAD_THR[it], PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                           "accepted bins;p_{T}^{reco,corr} [GeV/#it{c}];centrality-weighted counts",
                                           nbins_meas, bin_meas_edges);

          hRecoCorr_reject[it][ip] = new TH1D(Form("hRecoPtCorr_REJECTED_ptlead%.0f_%s_%s_%s", PTLEAD_THR[it], PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                             "rejected bins;p_{T}^{reco,corr} [GeV/#it{c}];centrality-weighted counts",
                                             nbins_meas, bin_meas_edges);
          hRecoCorr_entries[it][ip] = new TH1D(Form("hRecoPtCorr_ENTRIES_ptlead%.0f_%s_%s_%s", PTLEAD_THR[it], PtHatName(ip).Data(), rname.c_str(), cname.c_str()),
                                            "raw entries;p_{T}^{reco,corr} [GeV/#it{c}];raw entries",
                                            nbins_meas, bin_meas_edges);                                   
        }
      }

      for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        const int ip = FindPtHatBin((double)xsecWeight);
        if (ip < 0) continue;

        const bool haveMC   = (mc_pt > 0.0f);
        const bool haveReco = (reco_pt > (Float_t)RECO_DUMMY_CUT);
        if (!haveMC) continue;
        if (!haveReco) continue;
        if (!(deltaR > 0.0f && deltaR < (Float_t)dRmax)) continue;
        if (!(reco_area >= (Float_t)areaMin)) continue;
        if (!(reco_neutral_fraction <= (Float_t)CUT_NEUTRAL_FRACTION)) continue;
        if (!(reco_trigger_match == kTRUE)) continue;

        const double wCent = (double)centralityWeight;

        // Inclusive MC and reco plots after the same matched/reco-quality/trigger selection.
        hMc_pthat[ip]->Fill((double)mc_pt, wCent);
        hReco_pthat[ip]->Fill((double)reco_pt, wCent);

        for (int it = 0; it < N_LEAD; ++it) {
          if (reco_pt_lead >= (Float_t)PTLEAD_THR[it]) {
            hRecoCorr_overlay[it][ip]->Fill((double)reco_pt_corr, wCent);
            hRecoCorr_pthat[it][ip]->Fill((double)reco_pt_corr, wCent);
            hRecoCorr_entries[it][ip]->Fill((double)reco_pt_corr);
          }
        }
      }

      // Build accepted/rejected versions. The decision is made per pThat bin and per reco_corr bin.
      // Same criterion as your current code: content/error > sqrt(10), using the unscaled histogram.
      for (int ip = 0; ip < kNPthatBins; ++ip) {
        for (int it = 0; it < N_LEAD; ++it) {
          for (int ix = 1; ix <= nbins_meas; ++ix) {
            const double c = hRecoCorr_pthat[it][ip]->GetBinContent(ix);
            const double e = hRecoCorr_pthat[it][ip]->GetBinError(ix);
            const double signif = (e > 0.0 ? c / e : 0.0);
            const bool keep = (signif > kMinSignif);

            if (keep) {
              hRecoCorr_keep[it][ip]->SetBinContent(ix, c);
              hRecoCorr_keep[it][ip]->SetBinError(ix, e);
            } else {
              hRecoCorr_reject[it][ip]->SetBinContent(ix, c);
              hRecoCorr_reject[it][ip]->SetBinError(ix, e);
            }
          }
        }
      }

      outPthat->cd();
      for (int ip = 0; ip < kNPthatBins; ++ip) {
        hMc_pthat[ip]->Write();
        hReco_pthat[ip]->Write();
        for (int it = 0; it < N_LEAD; ++it) {
          hRecoCorr_overlay[it][ip]->Write();
          hRecoCorr_pthat[it][ip]->Write();
          hRecoCorr_entries[it][ip]->Write();
          hRecoCorr_keep[it][ip]->Write();
          hRecoCorr_reject[it][ip]->Write();
        }
      }

      TString centLabel = cname.c_str();
      centLabel.ReplaceAll("CENT_", "");
      centLabel.ReplaceAll("MID_",  "");
      centLabel.ReplaceAll("PERI_", "");
      centLabel.ReplaceAll("_", "-");
      centLabel += " %";

      // 1) one overlay plot for MC pT, one for reco pT.
      TString pdfName;
      pdfName.Form("%s/pthat_overlay_mc_%s_%s.pdf", pdfCDir.Data(), rname.c_str(), cname.c_str());
      DrawOverlayPthat(hMc_pthat, "MC #it{p}_{T} spectra by #hat{p}_{T} bin", "#it{p}_{T}^{MC} (GeV/#it{c})", pdfName, R, centLabel, 0, 60);

      pdfName.Form("%s/pthat_overlay_reco_%s_%s.pdf", pdfCDir.Data(), rname.c_str(), cname.c_str());
      DrawOverlayPthat(hReco_pthat, "Reco #it{p}_{T} spectra by #hat{p}_{T} bin", "#it{p}_{T}^{reco} (GeV/#it{c})", pdfName, R, centLabel, 0, 60);

      // 2) reco_corr overlay for every pTlead threshold.
      for (int it = 0; it < N_LEAD; ++it) {
        TH1D* tmp[kNPthatBins];
        for (int ip = 0; ip < kNPthatBins; ++ip) tmp[ip] = hRecoCorr_overlay[it][ip];
        pdfName.Form("%s/pthat_overlay_recoCorr_ptlead%.0f_%s_%s.pdf", pdfCDir.Data(), PTLEAD_THR[it], rname.c_str(), cname.c_str());
        DrawOverlayPthat(tmp,
                         Form("Reco corrected #it{p}_{T} spectra by #hat{p}_{T} bin, #it{p}_{T}^{lead} #geq %.0f", PTLEAD_THR[it]),
                         "#it{p}_{T,reco}^{corr} (GeV/#it{c})",
                         pdfName, R, centLabel, 0, 60);
      }

      // 3) significance mask plots: one per pThat bin and pTlead threshold.
      for (int it = 0; it < N_LEAD; ++it) {
        TString maskDir; maskDir.Form("%s/significance_ptlead%.0f", pdfCDir.Data(), PTLEAD_THR[it]);
        gSystem->mkdir(maskDir.Data(), kTRUE);

        for (int ip = 0; ip < kNPthatBins; ++ip) {
          pdfName.Form("%s/significance_%s_ptlead%.0f_%s_%s.pdf",
                       maskDir.Data(), PtHatName(ip).Data(), PTLEAD_THR[it], rname.c_str(), cname.c_str());
        DrawSignificanceMaskPlot(hRecoCorr_pthat[it][ip],
                                hRecoCorr_keep[it][ip],
                                hRecoCorr_reject[it][ip],
                                hRecoCorr_entries[it][ip],
                                ip, pdfName, R, centLabel, PTLEAD_THR[it]);
        }
      }

      for (int ip = 0; ip < kNPthatBins; ++ip) {
        delete hMc_pthat[ip];
        delete hReco_pthat[ip];
        for (int it = 0; it < N_LEAD; ++it) {
          delete hRecoCorr_pthat[it][ip];
          delete hRecoCorr_keep[it][ip];
          delete hRecoCorr_reject[it][ip];
          delete hRecoCorr_overlay[it][ip];
          delete hRecoCorr_entries[it][ip];
        }
      }
    }
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. ROOT output: " << outfile << "\n";
  std::cout << "PDFs written under: " << pdfTop << "/R*/CENT_* or MID_* or PERI_*\n";
}
