// unfold_data.C
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

// Same iterations as closure
static const int kBayesIterDefault = 4;

// measured & truth binning (must be IDENTICAL to what you used for response)
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

// same centrality & R labels as before
static const vector<string> kCentralities = {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii        = {"R0.2", "R0.3", "R0.4"};

// pt_lead cuts (same as in closure)
static const double kPtLeadCuts[] = {0.0, 5.0, 7.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// ------------------ jet-quality cuts -------------------------

static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

/**
 * dataFile      : real data jets (with JetTree)
 * responseDir   : directory where you stored unfold_response_*.root
 * outDir        : where unfolded data spectra will be written
 * nIter         : Bayes iterations (default 4)
 */
void unfold_data(const char* dataFile,
                 const char* responseDir,
                 const char* outDir,
                 int nIter = kBayesIterDefault)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    cout << "[error] Cannot open data file: " << dataFile << endl;
    return;
  }

  cout << "\n=== Unfolding real data ===\n";
  cout << "  Data file     : " << dataFile << "\n";
  cout << "  Response dir  : " << responseDir << "\n";
  cout << "  Output dir    : " << outDir << "\n";
  cout << "  Bayes iters   : " << nIter << "\n\n";

  // loop R, centrality, ptlead
  for (const auto& R : kRadii) {
    // choose area cut from R
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02;
    else if (Rval < 0.35) areaMin = CUT_AREA_03;
    else                  areaMin = CUT_AREA_04;

    for (const auto& C : kCentralities) {

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fData->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree in data: " << treePath << " (skip)\n";
        continue;
      }

      cout << "\n>>> Tag: " << R << "  " << C << "\n";

      // enable only needed branches (RECO only!)
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("xsecWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      float reco_pt_corr=0, reco_pt_lead=0;
      bool  reco_trigger_match=false;
      float centralityWeight=1.0f, xsecWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("xsecWeight", &xsecWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "  [note] empty tree in data, skip.\n";
        continue;
      }

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag      = R + "_" + C + Form("_ptlead%.0f", cut);
        const string tagfile  = tag; // same as used when saving responses

        cout << "  - pTlead cut >= " << cut << " GeV (tag " << tag << ")\n";

        // ---------------- measured spectrum in DATA ----------------
        TH1D* hMeasData = new TH1D(
            ("hMeasData_"+tag).c_str(),
            ";p_{T}^{reco,corr} [GeV];counts",
            nbins_meas, bin_meas_edges);
        hMeasData->SetDirectory(0);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "    filling data: " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          const double w = (double)centralityWeight * (double)xsecWeight;

          if (!reco_trigger_match)               continue;
          if (reco_area < areaMin)               continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;
          if (reco_pt_lead < cut)                continue;

          hMeasData->Fill(reco_pt_corr, w);
        }
        cout << "\n    data integral = " << hMeasData->Integral(0,-1) << endl;

        // ---------------- load response & prior ----------------
        const string respPath =
          string(responseDir) + "/unfold_response_" + tagfile + ".root";
        TFile* fResp = TFile::Open(respPath.c_str(), "READ");
        if (!fResp || fResp->IsZombie()) {
          cout << "    [ERROR] cannot open response file: " << respPath
               << "  --> skipping this tag\n";
          delete hMeasData;
          if (fResp) { fResp->Close(); delete fResp; }
          continue;
        }

        // response object name must match what you saved in closure macro
        const string respName = "response_" + tag;
        RooUnfoldResponse* response =
          dynamic_cast<RooUnfoldResponse*>(fResp->Get(respName.c_str()));
        if (!response) {
          cout << "    [ERROR] cannot find RooUnfoldResponse '" << respName
               << "' in " << respPath << " --> skipping\n";
          delete hMeasData;
          fResp->Close(); delete fResp;
          continue;
        }

        // prior: reuse same explicit prior as in closure (for consistency)
        const string priorName = "hPrior_" + tag;
        TH1D* hPrior = dynamic_cast<TH1D*>(fResp->Get(priorName.c_str()));
        if (!hPrior) {
          cout << "    [WARN] prior '" << priorName
               << "' not found, using internal response prior instead.\n";
        }

        // ---------------- run unfolding ----------------
        cout << "    running RooUnfoldBayes with " << nIter << " iterations\n";
        RooUnfoldBayes unfold(response, hMeasData, nIter,
                              /*smoothit=*/false,
                              hPrior); // may be null, RooUnfold handles it

        TH1D* hUnfoldedRecoBins = (TH1D*)unfold.Hunfold(); // reco-binning
        hUnfoldedRecoBins->SetName(Form("UnfoldedRecoBins_%s", tag.c_str()));
        hUnfoldedRecoBins->SetDirectory(0);

        // Rebin to truth binning (this is what you'll use later)
        TH1D* hUnfoldedTruthBins = (TH1D*)hUnfoldedRecoBins->Rebin(
            nbins_truth,
            Form("UnfoldedTruthBins_%s", tag.c_str()),
            bin_truth_edges);
        hUnfoldedTruthBins->SetDirectory(0);

        cout << "    unfolded integral (truth bins) = "
             << hUnfoldedTruthBins->Integral(0,-1) << endl;

        // ---------------- save to output file ----------------
        const string outPath =
          string(outDir) + "/unfolded_data_" + tagfile + ".root";
        TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
          cout << "    [ERROR] cannot create output file: " << outPath << endl;
        } else {
          hMeasData->Write("hMeasData");
          hUnfoldedRecoBins->Write("hUnfoldedRecoBins");
          hUnfoldedTruthBins->Write("hUnfoldedTruthBins");

          // optional: save covariance
          TMatrixD cov = unfold.Eunfold(RooUnfold::kCovariance);
          cov.Write("covariance");

          fout->Close();
          cout << "    wrote unfolded spectra to: " << outPath << endl;
        }

        // cleanup
        delete hMeasData;
        delete hUnfoldedRecoBins;
        delete hUnfoldedTruthBins;
        fResp->Close(); delete fResp;
      } // ptlead cuts
    } // centralities
  } // radii

  fData->Close();
  delete fData;

  cout << "\nAll done. Unfolded data written to: " << outDir << endl;
}
