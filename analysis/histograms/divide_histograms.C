#include <TFile.h>
#include <TH1.h>
#include <TDirectory.h>
#include <iostream>

void divide_histograms() {

    // --- Input files ---
    TFile *f1 = TFile::Open("hists_old.root", "READ");
    TFile *f2 = TFile::Open("hists_new.root", "READ");

   

    // --- Setup ---
    std::string radius[]     = {"R0.2", "R0.3", "R0.4"};
    std::string centrality[] = {"CENT_0_10", "MID_20_40", "PERI_60_80"};

    int N  = sizeof(radius)     / sizeof(radius[0]);
    int N2 = sizeof(centrality) / sizeof(centrality[0]);

    // --- Output file ---
    TFile *fout = new TFile("hists_ratio.root", "RECREATE");

    // --- Loop over radius ---
    for (int i = 0; i < N; i++) {

        std::cout << "Processing radius: " << radius[i] << std::endl;

        TDirectory *d1 = (TDirectory*)f1->Get(radius[i].c_str());
        TDirectory *d2 = (TDirectory*)f2->Get(radius[i].c_str());

        // Creating same directory structure in output file
        fout->mkdir(radius[i].c_str());
        fout->cd(radius[i].c_str());

        // --- Loop over centrality ---
        for (int j = 0; j < N2; j++) {

            std::cout << "  Centrality: " << centrality[j] << std::endl;

            TDirectory *sd1 = (TDirectory*)d1->Get(centrality[j].c_str());
            TDirectory *sd2 = (TDirectory*)d2->Get(centrality[j].c_str());

           

            // Histogram name
            std::string hname = "h_reco_ptcorr_" + radius[i] + "_" + centrality[j];

            TH1 *h1 = (TH1*)sd1->Get(hname.c_str());
            TH1 *h2 = (TH1*)sd2->Get(hname.c_str());

            if (!h1 || !h2) {
                std::cerr << "Error: histogram " << hname << " not found!" << std::endl;
                continue;
            }

            // Make ratio histogram
            std::string ratio_name = "h_ratio_" + radius[i] + "_" + centrality[j];
            TH1 *h_ratio = (TH1*)h1->Clone(ratio_name.c_str());
            h_ratio->Divide(h1, h2);

            // Create centrality subdir
            gDirectory->mkdir(centrality[j].c_str());
            gDirectory->cd(centrality[j].c_str());

            // Save histogram
            h_ratio->Write();

            // Return to radius directory for next loop
            fout->cd(radius[i].c_str());

            std::cout << "    Saved: " << ratio_name << std::endl;
        } // loop over centrality
    } // loop over radius

    fout->Close();
    std::cout << "All done!" << std::endl;
}
