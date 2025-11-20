#include <iostream>
#include <vector>
#include <iomanip> // For setting output precision

// ROOT includes
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TInterpreter.h"

// DD4hep includes
#include "DD4hep/Objects.h"
#include "DDG4/Geant4Data.h"

void readHits_Full() {
    // ============================================================
    // 1. LOAD LIBRARIES AND GENERATE DICTIONARY
    // ============================================================
    
    // Load DD4hep libraries (required for I/O and class definitions)
    if (gSystem->Load("libDDCore") < 0 && gSystem->Load("libDD4hep") < 0) {
        std::cerr << "Error: Could not load DD4hep core library." << std::endl;
        return;
    }
    gSystem->Load("libDDG4");
    gSystem->Load("libDDG4IO"); // Crucial for StreamerInfo/Dictionaries

    // Generate Dictionary for the hits vector (required for TTree::SetBranchAddress)
    gInterpreter->GenerateDictionary("vector<dd4hep::sim::Geant4Calorimeter::Hit*>", 
                                     "vector;DD4hep/Objects.h;DDG4/Geant4Data.h");

    // ============================================================
    // 2. OPEN FILE AND GET TREE
    // ============================================================
    TFile* file = new TFile("testSHiPCalo.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file 'testSHiPCalo.root'!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("EVENT");
    if (!tree) {
        std::cerr << "Error: TTree 'EVENT' not found!" << std::endl;
        file->Close();
        return;
    }

    // ============================================================
    // 3. SETUP BRANCH ADDRESS
    // ============================================================
    // Define the pointer to the vector of hits. Initialize to nullptr.
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* hits = nullptr;
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* Hhits = nullptr;

    if (!tree->GetBranch("SplitCalHits")) {
         std::cerr << "Error: Branch 'SplitCalHits' not found in tree!" << std::endl;
         file->Close();
         return;
    }

    tree->SetBranchAddress("SplitCalHits", &hits);
    tree->SetBranchAddress("SHiPHCALHits", &Hhits);

    // ============================================================
    // 4. LOOP OVER EVENTS AND HITS
    // ============================================================
    int nEvents = tree->GetEntries();
    std::cout << "--- Starting analysis of " << nEvents << " events ---" << std::endl;
    std::cout << std::fixed << std::setprecision(5); // Set precision for coordinates/energy

    TH1F *h_nrj = new TH1F("h_nrj","h_nrj;Energy loss [MeV];#",100,0,100);
    TH1F *h_x = new TH1F("h_x","h_x;X position [mm];#",100,-1000,1000);
    TH1F *h_y = new TH1F("h_y","h_y;Y position [mm];#",100,-1000,1000);
    TH1F *h_z = new TH1F("h_z","h_z;Z position [mm];#",2000,-1000,2500);
    TH2F *h_xz = new TH2F("h_xz","h_xz;X [mm];Z [mm]",300,-1000,1000,300,-1000,2500);
    TH2F *h_yz = new TH2F("h_yz","h_yz;Y [mm];Z [mm]",100,0,100,100,-1000,1000);

    for (int i = 0; i < nEvents; ++i) {
        tree->GetEntry(i);

        // Safety check for null pointer (shouldn't happen if dictionary loaded correctly)
        if (hits == nullptr) {
            std::cerr << "Fatal Error: Hits vector is null in Event " << i << ". Aborting." << std::endl;
            break; 
        }

        int nHits = hits->size();
        int nHHits = Hhits->size();
       
        for (size_t j = 0; j < Hhits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = Hhits->at(j);
            double energy = hit->energyDeposit;
            double x = hit->position.x();
            double y = hit->position.y();
            double z = hit->position.z();
		
	    cout << "Z "<< z << endl;

	    h_nrj->Fill(hit->energyDeposit);
	    h_x->Fill(hit->position.x());
	    h_y->Fill(hit->position.y());
	    h_z->Fill(hit->position.z());
	    h_xz->Fill(hit->position.x(),hit->position.z());
	    h_yz->Fill(hit->position.y(),hit->position.z());
	}

        // Loop over hits in this event
        for (size_t j = 0; j < hits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = hits->at(j);

            // Extract necessary information
            double energy = hit->energyDeposit;
            double x = hit->position.x();
            double y = hit->position.y();
            double z = hit->position.z();
		
	    h_nrj->Fill(hit->energyDeposit);
	    h_x->Fill(hit->position.x());
	    h_y->Fill(hit->position.y());
	    h_z->Fill(hit->position.z());
	    h_xz->Fill(hit->position.x(),hit->position.z());
	    h_yz->Fill(hit->position.y(),hit->position.z());

        }
    }

    TCanvas *c_nrj = new TCanvas("c_nrj","c_nrj",800,600);
    h_nrj->Draw();    
    TCanvas *c_x = new TCanvas("c_x","c_x",800,600);	
    h_x->Draw();    
    TCanvas *c_y = new TCanvas("c_y","c_y",800,600);	
    h_y->Draw();    
    TCanvas *c_z = new TCanvas("c_z","c_z",800,600);	
    h_z->Draw();    
    TCanvas *c_xz = new TCanvas("c_xz","c_xz",800,600);	
    h_xz->Draw("COLZ");    
    TCanvas *c_yz = new TCanvas("c_yz","c_yz",800,600);	
    h_yz->Draw("COLZ");    


//    file->Close();
}
