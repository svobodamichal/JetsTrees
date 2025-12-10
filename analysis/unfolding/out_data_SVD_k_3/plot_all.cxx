#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TROOT.h>
#include <iostream>

void plot_all() {
    // Disable ROOT's automatic object management
    TH1::AddDirectory(kFALSE);
    
    // Open the ROOT file
    TFile *file = TFile::Open("unfolded_data.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open unfolded_data.root" << std::endl;
        return;
    }
    
    std::cout << "File opened successfully: unfolded_data.root" << std::endl;
    
    // Create output directory
    gSystem->mkdir("plots", kTRUE);
    std::cout << "Output directory: plots/" << std::endl;
    
    // Get list of all keys (directories) in the file
    TList *keyList = file->GetListOfKeys();
    if (!keyList) {
        std::cerr << "ERROR: Cannot get list of keys from file" << std::endl;
        return;
    }
    
    int nProcessed = 0;
    int nSkipped = 0;
    
    // Iterate through all keys
    TIter nextKey(keyList);
    TKey *key;
    
    while ((key = (TKey*)nextKey())) {
        // Only process directory objects
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl || !cl->InheritsFrom("TDirectory")) {
            continue;
        }
        
        TString dirName = key->GetName();
        std::cout << "\n=== Processing: " << dirName << " ===" << std::endl;
        
        // Navigate into the directory
        TDirectory *currentDir = file->GetDirectory(dirName);
        
        if (!currentDir) {
            std::cerr << "WARNING: Cannot access directory " << dirName << std::endl;
            nSkipped++;
            continue;
        }
        
        // Get histograms and immediately clone them to detach from file
        TH1 *hUnfolded_orig = (TH1*)currentDir->Get("hUnfoldedTruthBins");
        TH1 *hMeasured_orig = (TH1*)currentDir->Get("hMeasData");
        TH1 *hBackfolded_orig = (TH1*)currentDir->Get("hBackfolded");
        
        // Check if all histograms exist
        if (!hUnfolded_orig || !hMeasured_orig || !hBackfolded_orig) {
            if (!hUnfolded_orig) std::cerr << "  WARNING: Missing hUnfoldedTruthBins" << std::endl;
            if (!hMeasured_orig) std::cerr << "  WARNING: Missing hMeasData" << std::endl;
            if (!hBackfolded_orig) std::cerr << "  WARNING: Missing hBackfolded" << std::endl;
            std::cerr << "  Skipping " << dirName << std::endl;
            nSkipped++;
            continue;
        }
        
        // Clone immediately to work with copies
        TH1 *hUnfolded = (TH1*)hUnfolded_orig->Clone(Form("unf_%d", nProcessed));
        TH1 *hMeasured = (TH1*)hMeasured_orig->Clone(Form("meas_%d", nProcessed));
        TH1 *hBackfolded = (TH1*)hBackfolded_orig->Clone(Form("back_%d", nProcessed));
        
        std::cout << "  Found all required histograms" << std::endl;
        
        // Create canvas
        TCanvas *canvas = new TCanvas(Form("c_%d", nProcessed), 
                                     Form("Unfolding: %s", dirName.Data()), 
                                     1200, 1200);
        
        // Divide canvas into two pads
        TPad *pad1 = new TPad(Form("pad1_%d", nProcessed), "Main", 0.0, 0.30, 1.0, 1.0);
        TPad *pad2 = new TPad(Form("pad2_%d", nProcessed), "Ratio", 0.0, 0.0, 1.0, 0.30);
        
        pad1->SetBottomMargin(0.1);
        pad1->SetTopMargin(0.055);
        pad1->SetLeftMargin(0.12);
        pad1->SetRightMargin(0.05);
        
        
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.25);
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.05);
        
        pad1->Draw();
        pad2->Draw();
        
        // === TOP PAD ===
        pad1->cd();
        pad1->SetLogy(1);
        
        // Style
        hMeasured->SetMarkerStyle(20);
        hMeasured->SetMarkerColor(kBlack);
        hMeasured->SetLineColor(kBlack);
        hMeasured->SetMarkerSize(1.0);
        hMeasured->SetTitle(dirName.Data());
        hMeasured->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}^{reco,corr} [(GeV/c)^{-1}]");
        hMeasured->GetYaxis()->SetTitleSize(0.05);
        hMeasured->GetYaxis()->SetLabelSize(0.04);
        hMeasured->GetYaxis()->SetTitleOffset(1.1);
        hMeasured->GetXaxis()->SetTitle("#it{p}_{T}^{reco,corr} [GeV/c]");
        hMeasured->GetXaxis()->SetRangeUser(-20, 65);
        
        hBackfolded->SetMarkerStyle(24);
        hBackfolded->SetMarkerColor(kRed+1);
        hBackfolded->SetLineColor(kRed+1);
        hBackfolded->SetMarkerSize(1.0);
        
        hUnfolded->SetMarkerStyle(21);
        hUnfolded->SetMarkerColor(kBlue+1);
        hUnfolded->SetLineColor(kBlue+1);
        hUnfolded->SetMarkerSize(1.0);
        
        // Draw
        hMeasured->Draw("EP");
        hBackfolded->Draw("SAME EP");
        hUnfolded->Draw("SAME EP");
        
        // Legend
        //TLegend *legend = new TLegend(0.15, 0.70, 0.45, 0.88);
        TLegend *legend = new TLegend(0.65, 0.70, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->AddEntry(hMeasured, "Measured Data", "lep");
        legend->AddEntry(hBackfolded, "Backfolded", "lep");
        legend->AddEntry(hUnfolded, "Unfolded", "lep");
        legend->Draw();
        
        // === BOTTOM PAD ===
        pad2->cd();
        
        // Create ratio
        TH1 *hRatio = (TH1*)hBackfolded->Clone(Form("ratio_%d", nProcessed));
        hRatio->Divide(hMeasured);
        
        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.8);
        hRatio->SetMarkerColor(kBlack);
        hRatio->SetLineColor(kBlack);
        hRatio->SetTitle("");
        
        hRatio->GetYaxis()->SetTitle("Backf/Meas");
        hRatio->GetYaxis()->SetTitleSize(0.1);
        hRatio->GetYaxis()->SetTitleOffset(0.45);
        hRatio->GetYaxis()->SetLabelSize(0.10);
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->SetMinimum(0.5);
        hRatio->SetMaximum(1.5);
        hRatio->GetXaxis()->SetRangeUser(-20, 65);
        
        hRatio->GetXaxis()->SetTitle("#it{p}_{T}^{reco,corr} [GeV/c]");
        hRatio->GetXaxis()->SetTitleSize(0.08);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        hRatio->GetXaxis()->SetTitleOffset(1.1);
        
        hRatio->Draw("EP");
        
        // Reference line
        TLine *line = new TLine(-20, 1.0, 
                                70, 1.0);
        line->SetLineStyle(2);
        line->SetLineColor(kGray+2);
        line->Draw();
        
        // Save
        canvas->Update();
        TString outputPath = Form("plots/%s.pdf", dirName.Data());
        canvas->SaveAs(outputPath);
        std::cout << "  Saved: " << outputPath << std::endl;
        
        // Clean up canvas (handles all drawn objects)
        delete canvas;
        
        // Clean up our cloned histograms
        delete hUnfolded;
        delete hMeasured;
        delete hBackfolded;
        
        nProcessed++;
    }
    
    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Successfully processed: " << nProcessed << " directories" << std::endl;
    std::cout << "Skipped: " << nSkipped << " directories" << std::endl;
    std::cout << "\nDone! All plots saved to plots/ directory." << std::endl;
    
    // DON'T call file->Close() - just let it go out of scope
    // The file will be closed automatically when ROOT exits
    // Calling Close() causes crashes due to directory cleanup issues
}