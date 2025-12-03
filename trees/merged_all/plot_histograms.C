#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TObjArray.h"

void plot_histograms() {
    // List of ROOT files
    const char* files[] = {
        "embedding_pt3_5.root", "embedding_pt5_7.root", "embedding_pt7_9.root", 
        "embedding_pt9_11.root", "embedding_pt11_15.root", "embedding_pt15_20.root", 
        "embedding_pt20_25.root", "embedding_pt25_30.root", "embedding_pt30_40.root", 
        "embedding_pt40_50.root", "embedding_pt50_-1.root"
    };
    const int nFiles = 11;
    
    const char* branchName = "reco_pt";
    
    // Set better plotting style
    gStyle->SetOptStat(0);
    
    // Create canvas
    TCanvas* c = new TCanvas("c", "Comparison", 2000, 1000);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetLogy(); // Optional: log scale for better visibility
    
    // Legend
    TLegend* legend = new TLegend(0.70, 0.45, 0.90, 0.90);
    legend->SetTextSize(0.03);
    //legend->SetBorderSize(0);
    
    // Color palette
    int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, 
                    kCyan+1, kViolet, kTeal, kPink, kAzure, kYellow+2};
    
    // Store histograms using TObjArray (ROOT 5 compatible)
    TObjArray* histograms = new TObjArray();
    histograms->SetOwner(kTRUE);
    
    float yMax = 0;
    int nSuccess = 0;
    
    for (int iFile = 0; iFile < nFiles; ++iFile) {
        TFile* f = TFile::Open(files[iFile], "READ");
        if (!f || f->IsZombie()) {
            printf("Cannot open %s\n", files[iFile]);
            if (f) delete f;
            continue;
        }
        
        TDirectory* dir = (TDirectory*)f->Get("R0.3");
        if (!dir) {
            printf("Directory R0.3 not found in %s\n", files[iFile]);
            f->Close();
            delete f;
            continue;
        }
        
        TDirectory* subdir = (TDirectory*)dir->Get("PERI_60_80");
        if (!subdir) {
            printf("Directory MID_40_60 not found in %s\n", files[iFile]);
            f->Close();
            delete f;
            continue;
        }
        
        TTree* tree = (TTree*)subdir->Get("JetTree");
        if (!tree) {
            printf("TTree JetTree not found in %s\n", files[iFile]);
            f->Close();
            delete f;
            continue;
        }
        
        // Create histogram for this file
        TH1D* h = new TH1D(Form("h_%d", iFile), ";reco_pt;Counts", 200, -20, 120);
        h->SetLineColor(colors[iFile % 11]);
        h->SetLineWidth(2);
        h->SetTitle("R0.3, PERI 60-80");
        h->SetDirectory(0); // Detach from file
        
        float value = 0;
        tree->SetBranchAddress(branchName, &value);
        
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++j) {
            tree->GetEntry(j);
            h->Fill(value);
        }
        
        // Track maximum for y-axis range
        if (h->GetMaximum() > yMax) {
            yMax = h->GetMaximum();
        }
        
        histograms->Add(h);
        
        // Extract pt range from filename for legend
        TString filename = files[iFile];
        filename.ReplaceAll("embedding_pt", "");
        filename.ReplaceAll(".root", "");
        legend->AddEntry(h, Form("pt %s", filename.Data()), "l");
        
        nSuccess++;
        
        f->Close();
        delete f;
    }
    
    // Draw all histograms
    for (int idx = 0; idx < histograms->GetEntries(); ++idx) {
        TH1D* h = (TH1D*)histograms->At(idx);
        if (idx == 0) {
            h->SetMaximum(yMax * 1.3);
            h->SetMinimum(0.1); // For log scale
            h->GetXaxis()->SetTitle("reco_pt");
            h->GetYaxis()->SetTitle("Counts");
            h->GetYaxis()->SetTitleOffset(1.4);
            h->Draw("HIST");
        } else {
            h->Draw("HIST SAME");
        }
    }
    
    legend->Draw();
    c->Update();
    
    // Optional: save to file
     c->SaveAs("reco_pt_comparison.pdf");
    
    printf("Successfully plotted %d histograms\n", nSuccess);
}