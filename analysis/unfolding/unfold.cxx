#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

static const int kBayesIters[] = {1, 4, 5, 7};
static const int kNBayesIters = sizeof(kBayesIters)/sizeof(kBayesIters[0]);

static const double kTestFrac = 0.50;     // 50/50 split
static const uint32_t kSeed    = 12345;   // deterministic split

static const double kPtLeadCuts[] = {0.0, 5.0, 7.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// measured & truth binning
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60
};
static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static const vector<string> kCentralities = {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii        = {"R0.2", "R0.3", "R0.4"};

// ------------------ jet-quality cuts -------------------------

// Area cuts per jet radius
static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

void unfold(const char* inputFile,
            const char* outDir)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);

  TFile* fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "[error] Cannot open input file: " << inputFile << endl;
    return;
  }

  // loop R, centrality, ptlead
  for (const auto& R : kRadii) {
    // parse numeric R to choose area cut
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02; // ~0.2
    else if (Rval < 0.35) areaMin = CUT_AREA_03; // ~0.3
    else                  areaMin = CUT_AREA_04; // ~0.4

    for (const auto& C : kCentralities) {

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fin->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree: " << treePath << " (skip)\n";
        continue;
      }

      // speed up I/O: enable only needed branches
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("mc_pt", 1);
      tr->SetBranchStatus("mc_pt_lead", 1);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("xsecWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      // branch addresses
      float mc_pt=0, mc_pt_lead=0;
      float reco_pt_corr=0, reco_pt_lead=0;
      bool  reco_trigger_match=false;
      float centralityWeight=1.0f, xsecWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("mc_pt", &mc_pt);
      tr->SetBranchAddress("mc_pt_lead", &mc_pt_lead);
      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("xsecWeight", &xsecWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "[note] empty tree: " << treePath << " (skip)\n";
        continue;
      }

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag = R + "_" + C + Form("_ptlead%.0f", cut);

        // --- histos: train/test for closure ---
        TH1D* hMeasTrain = new TH1D(("hMeas_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueTrain = new TH1D(("hTrue_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        TH1D* hMeasTest  = (TH1D*)hMeasTrain->Clone(("hMeasTest_"+tag).c_str());
        TH1D* hTrueTest  = (TH1D*)hTrueTrain->Clone(("hTrueTest_"+tag).c_str());

        hMeasTrain->SetDirectory(0);
        hTrueTrain->SetDirectory(0);
        hMeasTest ->SetDirectory(0);
        hTrueTest ->SetDirectory(0);

        // --- response matrix ---
        TH2D* hRespRecoVsTruth = new TH2D(("hResp_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);
        hRespRecoVsTruth->SetDirectory(0);

        // --- prior (truth-only, using *all* jets passing truth-side cuts) ---
        TH1D* hPrior = new TH1D(("hPrior_"+tag).c_str(),
            ";mc p_{T} [GeV];prior",
            nbins_truth, bin_truth_edges);
        hPrior->SetDirectory(0);

        TRandom3 rng(kSeed);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "  ["<< tag <<"] " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          // weight (same for prior and response)
          const double w = (double)centralityWeight * (double)xsecWeight;

          // ----- fill prior: only MC-side cuts -----
          if (mc_pt > 0.0 && mc_pt_lead >= cut) {
            hPrior->Fill(mc_pt, w);
          }

          // ----- reco-side cuts for response & closure -----
          if (!reco_trigger_match) continue;
          if (reco_area < areaMin) continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;

          // dual ptlead cut (both reco & MC)
          if (!(reco_pt_lead >= cut && mc_pt_lead >= cut)) continue;

          const bool train = (rng.Uniform() > kTestFrac);
          if (train) {
            hRespRecoVsTruth->Fill(reco_pt_corr, mc_pt, w);
            hMeasTrain->Fill(reco_pt_corr, w);
            hTrueTrain->Fill(mc_pt,        w);
          } else {
            hMeasTest->Fill(reco_pt_corr, w);
            hTrueTest->Fill(mc_pt,        w);
          }
        }
        cout << endl;

cout << "==== Sanity check for tag " << tag << " ====" << endl;
cout << "  Integral(Truth train) = " << hTrueTrain->Integral(0, -1) << endl;
cout << "  Integral(Meas  train) = " << hMeasTrain->Integral(0, -1) << endl;
cout << "  Integral(Truth test ) = " << hTrueTest ->Integral(0, -1) << endl;
cout << "  Integral(Meas  test ) = " << hMeasTest ->Integral(0, -1) << endl;
cout << "  Integral(Prior      ) = " << hPrior    ->Integral(0, -1) << endl;



        double intPrior = hPrior->Integral();
        double intTrue  = hTrueTrain->Integral();
        if (intPrior > 0.0 && intTrue > 0.0) {
          hPrior->Scale(intTrue / intPrior);
}

        // --- build response ---
        RooUnfoldResponse response(hMeasTrain, hTrueTrain, hRespRecoVsTruth);
        response.SetName(("response_"+tag).c_str());

        // --- unfolding with explicit prior ---
        vector<TH1D*> unfolded(kNBayesIters, nullptr);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          // RooUnfoldBayes(res, meas, nIter, smoothit=false, prior)
          RooUnfoldBayes u(&response, hMeasTest,
                           kBayesIters[ib], false, hPrior);
          TH1D* hunf = (TH1D*)u.Hunfold();
          hunf->SetDirectory(nullptr);
          hunf->SetName(Form("Unfolded_%s_iter%d", tag.c_str(), kBayesIters[ib]));
          unfolded[ib] = hunf;
        }

        // Rebin unfolded to truth binning (for plotting & ratio)
        vector<TH1D*> unfoldedTruth(kNBayesIters, nullptr);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          if (!unfolded[ib]) continue;
          unfoldedTruth[ib] = (TH1D*)unfolded[ib]->Rebin(
              nbins_truth,
              Form("UnfTruthBins_%d_%s", kBayesIters[ib], tag.c_str()),
              bin_truth_edges);
          unfoldedTruth[ib]->SetDirectory(nullptr);
        }

        // ========================= PLOTTING =========================
        TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1200);
        c->Divide(1,2);

        // ---------- top pad: shapes in TRUTH binning ----------
        c->cd(1);
        gPad->SetLogy();

        TH1D* hT_plot = (TH1D*)hTrueTest->Clone(("hTrueW_"+tag).c_str());
      //  hT_plot->Scale(1.0, "width");
        hT_plot->SetMarkerStyle(20);
        hT_plot->SetLineColor(kBlack);

        TH1D* hM_truth = (TH1D*)hMeasTest->Rebin(
            nbins_truth,
            ("hMeasTruthBins_"+tag).c_str(),
            bin_truth_edges);
      //  hM_truth->Scale(1.0, "width");
        hM_truth->SetMarkerStyle(24);
        hM_truth->SetLineColor(kBlue+1);

        hT_plot->GetXaxis()->SetTitle("p_{T} [GeV]");
        hT_plot->Draw("E1");
        hM_truth->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.58,0.58,0.88,0.88);
        leg->AddEntry(hT_plot,  "Truth (test)",    "lp");
        leg->AddEntry(hM_truth, "Measured (test)", "lp");

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* w = unfoldedTruth[ib];
          if (!w) continue;
       //   w->Scale(1.0, "width");
          w->SetMarkerStyle(20);
          w->SetMarkerColor(kMagenta + ib);
          w->SetLineColor(kMagenta + ib);
          w->Draw("E1 SAME");
          leg->AddEntry(w, Form("Bayes %d it.", kBayesIters[ib]), "lp");
        }
        leg->Draw();

        // ---------- bottom pad: ratio unfolded / truth ----------
        c->cd(2);
        gPad->SetGridy();  // optional, just cosmetic

        TH1D* firstRatio = nullptr;

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* hunf_reb = unfoldedTruth[ib];
          if (!hunf_reb) continue;

          TH1D* r = (TH1D*)hunf_reb->Clone(
              Form("ratio_%d_%s", kBayesIters[ib], tag.c_str()));
          r->Divide(hTrueTest);       // bin–by–bin unfolded / truth

          if (!firstRatio) {
            firstRatio = r;
            firstRatio->SetTitle("");
            firstRatio->GetYaxis()->SetTitle("Unfolded / Truth");
            firstRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
            firstRatio->Draw("E1");   // defines axes
          } else {
            r->Draw("E1 SAME");
          }
        }

        // draw guide lines at 0.95 and 1.05 if we actually drew something
        if (firstRatio) {
          double xmin = firstRatio->GetXaxis()->GetXmin();
          double xmax = firstRatio->GetXaxis()->GetXmax();

          TLine* ln1 = new TLine(xmin, 1.05, xmax, 1.05);
          ln1->SetLineStyle(2); ln1->SetLineColor(kGray+1); ln1->Draw();

          TLine* ln2 = new TLine(xmin, 0.95, xmax, 0.95);
          ln2->SetLineStyle(2); ln2->SetLineColor(kGray+1); ln2->Draw();
        }

        // ---------- save ----------
        const string tagfile = R + "_" + C + Form("_ptlead%.0f", cut);
        const string rootPath = string(outDir) + "/unfold_response_" + tagfile + ".root";
        const string pdfPath  = string(outDir) + "/closure_" + tagfile + ".pdf";

        TFile* outf = TFile::Open(rootPath.c_str(), "RECREATE");
        hRespRecoVsTruth->Write();
        response.Write();
        hMeasTrain->Write(); hTrueTrain->Write();
        hMeasTest->Write();  hTrueTest->Write();
        hPrior->Write();     // save the prior used
        for (auto* u : unfolded) if (u) u->Write();
        outf->Close();

        c->SaveAs(pdfPath.c_str());

        // ---------- cleanup ----------
        delete c;
        delete hT_plot;
        delete hM_truth;
        delete hRespRecoVsTruth;
        delete hMeasTrain;
        delete hTrueTrain;
        delete hMeasTest;
        delete hTrueTest;
        delete hPrior;
        delete leg;
        for (auto* u  : unfolded)      delete u;
        for (auto* ut : unfoldedTruth) delete ut;
      } // ptlead cuts
    } // centralities
  } // radii

  cout << "Done. Outputs in: " << outDir << endl;
  fin->Close();
  delete fin;
}