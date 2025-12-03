#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TCollection.h"

#include <iostream>
#include <string>
#include <cstdio>

// pTlead thresholds (GeV)
const int N_PLEAD = 4;
double PLEAD_CUTS[N_PLEAD] = {0.0, 5.0, 7.0, 9.0};

void make_efficienciesH(const char *infile  = "embedding_mergedH.root",
                       const char *outfile = "efficienciesH.root")
{
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

  // ========== loop over R directories (R0.2, R0.3, ...) ==========
  TIter nextR(fin->GetListOfKeys());
  TKey *keyR = 0;

  while ((keyR = (TKey*) nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.size() < 2 || rname[0] != 'R') continue;

    // parse R value from e.g. "R0.2"
    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) {
      std::cerr << "  Warning: cannot parse R from '" << rname
                << "'. Skipping." << std::endl;
      continue;
    }

    TDirectoryFile *dirR = (TDirectoryFile*) keyR->ReadObj();
    if (!dirR) continue;

    std::cout << "Radius directory: " << rname << " (R=" << R << ")" << std::endl;

    fout->cd();
    TDirectory *outR = fout->mkdir(rname.c_str());
    if (!outR) outR = fout->GetDirectory(rname.c_str());

    // ========== loop over centrality subdirectories ==========
    TIter nextC(dirR->GetListOfKeys());
    TKey *keyC = 0;

    while ((keyC = (TKey*) nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      TDirectoryFile *dirC = (TDirectoryFile*) keyC->ReadObj();
      if (!dirC) continue;

      std::cout << "  Centrality: " << cname << std::endl;

      // fetch histograms we need from the old output
      TH2D *hDen2          = (TH2D*) dirC->Get("den_ptcorr_vs_ptlead");
      TH2D *hNum2          = (TH2D*) dirC->Get("num_ptcorr_vs_ptlead");
      TH2D *hRecoMatched2  = (TH2D*) dirC->Get("reco_matched_ptcorr_vs_ptlead");
      TH1D *hMc1           = (TH1D*) dirC->Get("mc_pt");
      TH2D *hRecoMc2       = (TH2D*) dirC->Get("recoptcorr_vs_mcpt");

      if (!hDen2 || !hNum2 || !hRecoMatched2 || !hMc1 || !hRecoMc2) {
        std::cerr << "    Warning: missing histos in "
                  << rname << "/" << cname << "  -> skipping." << std::endl;
        continue;
      }

      outR->cd();
      TDirectory *outC = outR->mkdir(cname.c_str());
      if (!outC) outC = outR->GetDirectory(cname.c_str());
      outC->cd();

      // -------- bin info (for pTlead axis etc.) ----------
      int    ny_lead    = hDen2->GetYaxis()->GetNbins();
      double y_lead_min = hDen2->GetYaxis()->GetXmin();
      double y_lead_max = hDen2->GetYaxis()->GetXmax();

      // ========== loop over pTlead thresholds ==========
      for (int iL = 0; iL < N_PLEAD; ++iL) {
        double leadCut = PLEAD_CUTS[iL];

        // if cut is beyond axis max, skip
        if (leadCut >= y_lead_max) {
          std::cout << "    PLEAD cut " << leadCut
                    << " >= max axis " << y_lead_max
                    << " -> skipping." << std::endl;
          continue;
        }

        int leadInt = (int)(leadCut + 0.5);
        char leadDirName[32];
        sprintf(leadDirName, "PLEAD_%d", leadInt);

        outC->cd();
        TDirectory *outL = outC->mkdir(leadDirName);
        if (!outL) outL = outC->GetDirectory(leadDirName);
        outL->cd();

        std::cout << "    pTlead >= " << leadCut
                  << "  (dir: " << leadDirName << ")" << std::endl;

        // find y-bin corresponding to pTlead >= leadCut
        int ybin_min = hDen2->GetYaxis()->FindBin(leadCut + 1e-6);
        int ybin_max = ny_lead;

        //=============================
        // 1) Trigger efficiency
        //=============================
        // projections of Den and Num on reco pTcorr axis
        char nameDenTrig[256], nameNumTrig[256];
        sprintf(nameDenTrig, "den_trig_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);
        sprintf(nameNumTrig, "num_trig_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);

        TH1D *hDenTrig = hDen2->ProjectionX(nameDenTrig, ybin_min, ybin_max, "e");
        TH1D *hNumTrig = hNum2->ProjectionX(nameNumTrig, ybin_min, ybin_max, "e");

        // binomial efficiency: hTrigEff = Num / Den
        char nameTrigEff[256];
        sprintf(nameTrigEff, "trigEff_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);

        TH1D *hTrigEff = (TH1D*) hNumTrig->Clone(nameTrigEff);
        hTrigEff->Reset();
        hTrigEff->Divide(hNumTrig, hDenTrig, 1.0, 1.0, "b");
        hTrigEff->SetTitle("Trigger efficiency vs reco p_{T}^{corr}");

        //=============================
        // 2) Purity and fake rate
        //=============================
        // Denominator: all reco jets (same Den as trigger)
        // Numerator: matched reco jets (from reco_matched_ptcorr_vs_ptlead)
        char nameDenPur[256], nameNumPur[256];
        sprintf(nameDenPur, "den_purity_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);
        sprintf(nameNumPur, "num_purity_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);

        TH1D *hDenPur = hDen2->ProjectionX(nameDenPur, ybin_min, ybin_max, "e");
        TH1D *hNumPur = hRecoMatched2->ProjectionX(nameNumPur, ybin_min, ybin_max, "e");

        char namePurity[256], nameFake[256];
        sprintf(namePurity, "purity_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);
        sprintf(nameFake,   "fakeRate_ptcorr_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);

        TH1D *hPurity = (TH1D*) hNumPur->Clone(namePurity);
        hPurity->Reset();
        hPurity->Divide(hNumPur, hDenPur, 1.0, 1.0, "b");
        hPurity->SetTitle("Purity vs reco p_{T}^{corr} (matched / all reco)");

        // fake rate = 1 - purity
        TH1D *hFake = (TH1D*) hPurity->Clone(nameFake);
        hFake->SetTitle("Fake rate vs reco p_{T}^{corr}");
        int nb = hFake->GetNbinsX();
        for (int ib = 1; ib <= nb; ++ib) {
          double p = hPurity->GetBinContent(ib);
          double e = hPurity->GetBinError(ib);
          hFake->SetBinContent(ib, 1.0 - p);
          hFake->SetBinError(ib, e); // same absolute error
        }

        //=============================
        // 3) Matching efficiency vs MC pT
        //=============================
        // This does NOT depend on pTlead in this old-output scheme.
        // Denominator: hMc1 (all MC jets)
        // Numerator: projection of reco_mc on MC pt axis
        char nameMcMatch[256], nameMatchEff[256];
        sprintf(nameMcMatch,   "mc_matched_pt_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);
        sprintf(nameMatchEff,  "matchEff_mcpt_%s_%s_%s",
                rname.c_str(), cname.c_str(), leadDirName);

        // proj along X (MC pT) over all reco pTcorr bins
        int ybin_min_mc = 1;
        int ybin_max_mc = hRecoMc2->GetYaxis()->GetNbins();
        TH1D *hMcMatched = hRecoMc2->ProjectionX(nameMcMatch,
                                                 ybin_min_mc, ybin_max_mc, "e");

        TH1D *hMatchEff = (TH1D*) hMcMatched->Clone(nameMatchEff);
        hMatchEff->Reset();
        hMatchEff->Divide(hMcMatched, hMc1, 1.0, 1.0, "b");
        hMatchEff->SetTitle("Matching efficiency vs MC p_{T}");

        // (Do NOT call Write() here; we let fout->Write() at the end handle it.)
      } // end loop pTlead
    } // end loop centrality
  } // end loop radius

  fout->Write();   // single global write => no duplicate keys
  fout->Close();
  fin->Close();

  std::cout << "Done. Efficiencies written to " << outfile << std::endl;
}
