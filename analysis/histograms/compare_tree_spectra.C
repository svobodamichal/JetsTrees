#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <cmath>

const double CUT_AREA_02 = 0.07;
const double CUT_AREA_03 = 0.20;
const double CUT_AREA_04 = 0.40;
const double CUT_NEUTRAL_FRACTION = 0.95;
const double RECO_DUMMY_CUT = -500.0;

const int N_R = 3;
const double RVALS[N_R] = {0.2, 0.3, 0.4};

const int N_CENT = 3;
const char* CENTS[N_CENT] = {
  "CENT_0_10",
  "MID_20_40",
  "PERI_60_80"
};

const int N_LEAD = 4;
const double PTLEAD[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

const int nbins_meas = 24;
const double bin_meas_edges[nbins_meas + 1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,
  10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60
};

const int nbins_fine = 100;
const double xmin_fine = -40.0;
const double xmax_fine = 60.0;

double areaMinForR(double R)
{
  if (R < 0.25) return CUT_AREA_02;
  if (R < 0.35) return CUT_AREA_03;
  return CUT_AREA_04;
}

TH1D* makeSpectrumFromTree(const char* fileName,
                           double R,
                           const char* cent,
                           double ptlead,
                           const char* hname,
                           bool isEmbedding,
                           bool fineBinning)
{
  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << fileName << std::endl;
    return 0;
  }

  TString treePath;
  treePath.Form("R%.1f/%s/JetTree", R, cent);

  TTree* t = (TTree*)f->Get(treePath.Data());
  if (!t) {
    std::cerr << "Missing tree " << treePath << " in " << fileName << std::endl;
    f->Close();
    return 0;
  }

  Float_t centralityWeight = 1.0;
  Float_t xsecWeight = 1.0;

  Float_t reco_pt = -999.0;
  Float_t reco_pt_corr = -999.0;
  Float_t reco_area = -999.0;
  Float_t reco_neutral_fraction = 999.0;
  Float_t reco_pt_lead = -999.0;
  Bool_t reco_trigger_match = kFALSE;

  t->SetBranchAddress("centralityWeight", &centralityWeight);
  t->SetBranchAddress("reco_pt", &reco_pt);
  t->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
  t->SetBranchAddress("reco_area", &reco_area);
  t->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);
  t->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
  t->SetBranchAddress("reco_trigger_match", &reco_trigger_match);

  if (isEmbedding && t->GetBranch("xsecWeight")) {
    t->SetBranchAddress("xsecWeight", &xsecWeight);
  }

  TH1D* h = 0;
    if (fineBinning) {
      h = new TH1D(hname, hname, nbins_fine, xmin_fine, xmax_fine);
    } else {
      h = new TH1D(hname, hname, nbins_meas, bin_meas_edges);
    }
  h->Sumw2();

  const double areaMin = areaMinForR(R);
  const Long64_t nentries = t->GetEntries();

  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);

    if (!(reco_pt > RECO_DUMMY_CUT)) continue;
    if (!(reco_area >= areaMin)) continue;
    if (!(reco_neutral_fraction <= CUT_NEUTRAL_FRACTION)) continue;
    if (!(reco_trigger_match == kTRUE)) continue;
    if (!(reco_pt_lead >= ptlead)) continue;

    double w = centralityWeight;
    if (isEmbedding) w *= xsecWeight;

    h->Fill(reco_pt_corr, w);
  }

  h->SetDirectory(0);
  f->Close();

  // bin-width correction
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double bw = h->GetBinWidth(ib);
    if (bw <= 0.0) continue;
    h->SetBinContent(ib, h->GetBinContent(ib) / bw);
    h->SetBinError(ib, h->GetBinError(ib) / bw);
  }

  // integral normalization after bin-width correction
  const double integral = h->Integral("width");
  if (integral > 0.0) h->Scale(1.0 / integral);

  return h;
}

void compare_tree_spectra(const char* treeDir = "../../trees",
                          const char* outDir = "comparison_tree_spectra")
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  const char* dataFile = "data_merged.root";

  const int N_SET = 4;
  const char* embFiles[N_SET] = {
    "embedding_merged_MCReco1p5.root",
    "embedding_merged_MCReco2p0.root",
    "embedding_merged_NoRecoVeto.root",
    "embedding_merged_NoVeto.root"
  };

  const char* embLabels[N_SET] = {
    "MCReco1p5",
    "MCReco2p0",
    "NoRecoVeto",
    "NoVeto"
  };

  const int colors[N_SET] = {
    kRed + 1,
    kBlue + 1,
    kGreen + 2,
    kMagenta + 2
  };

  const int markers[N_SET] = {
    21, 22, 23, 33
  };

  TString dataPath;
  dataPath.Form("%s/%s", treeDir, dataFile);

  for (int ir = 0; ir < N_R; ++ir) {
    const double R = RVALS[ir];

    for (int ic = 0; ic < N_CENT; ++ic) {
      const char* cent = CENTS[ic];

      TString subDir;
      subDir.Form("%s/R%.1f/%s", outDir, R, cent);
      gSystem->mkdir(subDir.Data(), kTRUE);

      for (int il = 0; il < N_LEAD; ++il) {
        const double ptlead = PTLEAD[il];

        TH1D* hData = makeSpectrumFromTree(
          dataPath.Data(), R, cent, ptlead,
          Form("hData_R%.1f_%s_pTlead%.0f", R, cent, ptlead),
          false,
          false
        );

        if (!hData) continue;

        hData->SetMarkerStyle(20);
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);

        TH1D* hEmb[N_SET] = {0};
        TH1D* hRatio[N_SET] = {0};

        for (int is = 0; is < N_SET; ++is) {
          TString embPath;
          embPath.Form("%s/%s", treeDir, embFiles[is]);

          hEmb[is] = makeSpectrumFromTree(
            embPath.Data(), R, cent, ptlead,
            Form("hEmb_%s_R%.1f_%s_pTlead%.0f", embLabels[is], R, cent, ptlead),
            true,
            false
          );

          if (!hEmb[is]) continue;

          hEmb[is]->SetMarkerStyle(markers[is]);
          hEmb[is]->SetMarkerColor(colors[is]);
          hEmb[is]->SetLineColor(colors[is]);

          hRatio[is] = (TH1D*)hEmb[is]->Clone(
            Form("hRatio_%s_R%.1f_%s_pTlead%.0f", embLabels[is], R, cent, ptlead)
          );
          hRatio[is]->Divide(hData);
        }

        TCanvas* c = new TCanvas("c", "c", 850, 850);

        TPad* pTop = new TPad("pTop", "pTop", 0.0, 0.30, 1.0, 1.0);
        TPad* pBot = new TPad("pBot", "pBot", 0.0, 0.00, 1.0, 0.30);

        pTop->SetBottomMargin(0.02);
        pTop->SetLogy();

        pBot->SetTopMargin(0.03);
        pBot->SetBottomMargin(0.32);

        pTop->Draw();
        pBot->Draw();

        pTop->cd();

        hData->SetTitle("");
        hData->GetXaxis()->SetLabelSize(0);
        hData->GetXaxis()->SetTitleSize(0);
        hData->GetXaxis()->SetRangeUser(-10.0, 60.0);
        hData->GetYaxis()->SetTitle("normalized spectra");
        hData->GetYaxis()->SetTitleOffset(1.25);
        hData->SetMinimum(1e-6);
        hData->SetMaximum(1.0);
        hData->Draw("E1");

        for (int is = 0; is < N_SET; ++is) {
          if (hEmb[is]) hEmb[is]->Draw("E1 SAME");
        }

        TLegend* leg = new TLegend(0.56, 0.62, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(hData, "data", "lep");

        for (int is = 0; is < N_SET; ++is) {
          if (hEmb[is]) leg->AddEntry(hEmb[is], embLabels[is], "lep");
        }

        leg->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.040);

        TString centLabel = cent;
        centLabel.ReplaceAll("CENT_", "");
        centLabel.ReplaceAll("MID_", "");
        centLabel.ReplaceAll("PERI_", "");
        centLabel.ReplaceAll("_", "-");

        tex.DrawLatex(0.18, 0.23, Form("R = %.1f, %s %%", R, centLabel.Data()));
        tex.DrawLatex(0.18, 0.17, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", ptlead));
        tex.DrawLatex(0.18, 0.11, "Integral-normalized reco spectra");

        pBot->cd();

        TH1D* frame = 0;
        for (int is = 0; is < N_SET; ++is) {
          if (hRatio[is]) {
            frame = hRatio[is];
            break;
          }
        }

        if (frame) {
          frame->SetTitle("");
          frame->GetXaxis()->SetTitle("#it{p}_{T,reco}^{corr} (GeV/#it{c})");
          frame->GetYaxis()->SetTitle("embedding / data");
          frame->GetXaxis()->SetRangeUser(-10.0, 60.0);
          frame->GetYaxis()->SetRangeUser(0.0, 2.0);

          frame->GetXaxis()->SetTitleSize(0.11);
          frame->GetXaxis()->SetLabelSize(0.09);
          frame->GetYaxis()->SetTitleSize(0.09);
          frame->GetYaxis()->SetLabelSize(0.08);
          frame->GetYaxis()->SetTitleOffset(0.55);
          frame->GetYaxis()->SetNdivisions(505);

          frame->Draw("E1");

          for (int is = 0; is < N_SET; ++is) {
            if (hRatio[is] && hRatio[is] != frame) hRatio[is]->Draw("E1 SAME");
          }

          TLine* line = new TLine(-10.0, 1.0, 60.0, 1.0);
          line->SetLineStyle(2);
          line->Draw("SAME");
        }

        TString outPng;
        outPng.Form("%s/compare_treeSpectra_R%.1f_%s_pTlead%.0f.png",
                    subDir.Data(), R, cent, ptlead);

        c->SaveAs(outPng.Data());

        delete c;
        delete leg;
        delete hData;

        // =====================================================
        // Fine-binning version: 100 bins from -40 to 60 GeV
        // PNG only
        // =====================================================
        {
          TH1D* hDataFine = makeSpectrumFromTree(
            dataPath.Data(), R, cent, ptlead,
            Form("hDataFine_R%.1f_%s_pTlead%.0f", R, cent, ptlead),
            false,
            true
          );

          if (!hDataFine) continue;

          hDataFine->SetMarkerStyle(20);
          hDataFine->SetMarkerColor(kBlack);
          hDataFine->SetLineColor(kBlack);

          TH1D* hEmbFine[N_SET] = {0};
          TH1D* hRatioFine[N_SET] = {0};

          for (int is = 0; is < N_SET; ++is) {
            TString embPath;
            embPath.Form("%s/%s", treeDir, embFiles[is]);

            hEmbFine[is] = makeSpectrumFromTree(
              embPath.Data(), R, cent, ptlead,
              Form("hEmbFine_%s_R%.1f_%s_pTlead%.0f", embLabels[is], R, cent, ptlead),
              true,
              true
            );

            if (!hEmbFine[is]) continue;

            hEmbFine[is]->SetMarkerStyle(markers[is]);
            hEmbFine[is]->SetMarkerColor(colors[is]);
            hEmbFine[is]->SetLineColor(colors[is]);

            hRatioFine[is] = (TH1D*)hEmbFine[is]->Clone(
              Form("hRatioFine_%s_R%.1f_%s_pTlead%.0f", embLabels[is], R, cent, ptlead)
            );

            hRatioFine[is]->Divide(hDataFine);
          }

          TCanvas* cFine = new TCanvas("cFine", "cFine", 850, 850);

          TPad* pTopFine = new TPad("pTopFine", "pTopFine", 0.0, 0.30, 1.0, 1.0);
          TPad* pBotFine = new TPad("pBotFine", "pBotFine", 0.0, 0.00, 1.0, 0.30);

          pTopFine->SetBottomMargin(0.02);
          pTopFine->SetLogy();

          pBotFine->SetTopMargin(0.03);
          pBotFine->SetBottomMargin(0.32);

          pTopFine->Draw();
          pBotFine->Draw();

          pTopFine->cd();

          hDataFine->SetTitle("");
          hDataFine->GetXaxis()->SetLabelSize(0);
          hDataFine->GetXaxis()->SetTitleSize(0);
          hDataFine->GetXaxis()->SetRangeUser(-40.0, 60.0);
          hDataFine->GetYaxis()->SetTitle("normalized spectra");
          hDataFine->GetYaxis()->SetTitleOffset(1.25);
          hDataFine->SetMinimum(1e-6);
          hDataFine->SetMaximum(1.0);
          hDataFine->Draw("E1");

          for (int is = 0; is < N_SET; ++is) {
            if (hEmbFine[is]) hEmbFine[is]->Draw("E1 SAME");
          }

          TLegend* legFine = new TLegend(0.56, 0.62, 0.88, 0.88);
          legFine->SetBorderSize(0);
          legFine->SetFillStyle(0);
          legFine->SetTextSize(0.035);
          legFine->AddEntry(hDataFine, "data", "lep");

          for (int is = 0; is < N_SET; ++is) {
            if (hEmbFine[is]) legFine->AddEntry(hEmbFine[is], embLabels[is], "lep");
          }

          legFine->Draw();

          TLatex texFine;
          texFine.SetNDC();
          texFine.SetTextFont(42);
          texFine.SetTextSize(0.040);

          TString centLabelFine = cent;
          centLabelFine.ReplaceAll("CENT_", "");
          centLabelFine.ReplaceAll("MID_", "");
          centLabelFine.ReplaceAll("PERI_", "");
          centLabelFine.ReplaceAll("_", "-");

          texFine.DrawLatex(0.18, 0.23, Form("R = %.1f, %s %%", R, centLabelFine.Data()));
          texFine.DrawLatex(0.18, 0.17, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", ptlead));
          texFine.DrawLatex(0.18, 0.11, "Fine binning, integral-normalized");

          pBotFine->cd();

          TH1D* frameFine = 0;

          for (int is = 0; is < N_SET; ++is) {
            if (hRatioFine[is]) {
              frameFine = hRatioFine[is];
              break;
            }
          }

          if (frameFine) {
            frameFine->SetTitle("");
            frameFine->GetXaxis()->SetTitle("#it{p}_{T,reco}^{corr} (GeV/#it{c})");
            frameFine->GetYaxis()->SetTitle("embedding / data");
            frameFine->GetXaxis()->SetRangeUser(-40.0, 60.0);
            frameFine->GetYaxis()->SetRangeUser(0.0, 2.0);

            frameFine->GetXaxis()->SetTitleSize(0.11);
            frameFine->GetXaxis()->SetLabelSize(0.09);
            frameFine->GetYaxis()->SetTitleSize(0.09);
            frameFine->GetYaxis()->SetLabelSize(0.08);
            frameFine->GetYaxis()->SetTitleOffset(0.55);
            frameFine->GetYaxis()->SetNdivisions(505);

            frameFine->Draw("E1");

            for (int is = 0; is < N_SET; ++is) {
              if (hRatioFine[is] && hRatioFine[is] != frameFine) {
                hRatioFine[is]->Draw("E1 SAME");
              }
            }

            TLine* lineFine = new TLine(-40.0, 1.0, 60.0, 1.0);
            lineFine->SetLineStyle(2);
            lineFine->Draw("SAME");
          }

          TString outPngFine;
          outPngFine.Form("%s/compare_treeSpectra_FINE_R%.1f_%s_pTlead%.0f.png",
                          subDir.Data(), R, cent, ptlead);

          cFine->SaveAs(outPngFine.Data());

          delete cFine;
          delete legFine;
          delete hDataFine;

          for (int is = 0; is < N_SET; ++is) {
            delete hEmbFine[is];
            delete hRatioFine[is];
          }
        }

        for (int is = 0; is < N_SET; ++is) {
          delete hEmb[is];
          delete hRatio[is];
        }
      }
    }
  }

  std::cout << "Done. Output written to: " << outDir << std::endl;
}