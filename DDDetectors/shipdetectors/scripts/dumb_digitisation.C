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

void Reset_nhits(int *arr){

	for(int i=0;i<5;i++){
		arr[i] = 0;
	}

}

double vec_sum(const std::vector<double> v){
	double sum = 0;
	for(auto a:v){
		sum += a;
	}
	return sum;
}

double vec_av(const std::vector<double> v){
	return vec_sum(v)/static_cast<double>(v.size());
}

double vec_rms(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0; 
    }
    double sum_of_squares = 0.0;
    for (double value : data) {
        sum_of_squares += value * value;
    }
    double mean_of_squares = sum_of_squares / data.size();
    return TMath::Sqrt(mean_of_squares);
}


int DecodeLayer(double z){
	int layer = -9999;
	if(z > 0 && z < 100) layer = 0;
	if(z > 100 && z < 300) layer = 1;
	if(z > 300 && z < 400) layer = 2;
	if(z > 400 && z < 600) layer = 3;
	if(z > 600 && z < 800) layer = 4;
    return layer;
}



void setup_data_csv(std::string inparticle, int energy, int runid) {
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

    TString instr = "/eos/user/m/mclimesc/SPLITCAL/SplitCalPhysics/DD4HEP_PID/PID_NoSplitCal/"+inparticle+"_sample_%dGeV_%d.root";
    cout << Form(instr,energy) << endl; 
    TFile* file = new TFile(Form(instr,energy,runid), "READ");
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
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* Hhits = nullptr;

    tree->SetBranchAddress("SHiPHCALHits", &Hhits);

    // ============================================================
    // 4. LOOP OVER EVENTS AND HITS
    // ============================================================
    int nEvents = tree->GetEntries();
    std::cout << "--- Starting analysis of " << nEvents << " events ---" << std::endl;
    std::cout << std::fixed << std::setprecision(5); // Set precision for coordinates/energy

    std::stringstream filename;
    filename << "/eos/user/m/mclimesc/SPLITCAL/SplitCalPhysics/DD4HEP_PID/PID_NoSplitCal/" << inparticle << "_sample_" << energy << "GeV_"<< runid <<".csv";

    ofstream outfile;
    outfile.open(filename.str());
    
    outfile << "event,"; 

    for(int hcallayer=0;hcallayer<5;hcallayer++){

    	outfile << "avx_" << hcallayer << ","; 
    	outfile << "avy_" << hcallayer << ","; 
    	outfile << "rmsx_" << hcallayer << ","; 
    	outfile << "rmsy_" << hcallayer << ","; 
    	outfile << "nhits_" << hcallayer << ","; 
    	outfile << "sumenergydep_" << hcallayer << ","; 
    	outfile << "rmsenergydep_" << hcallayer << ","; 
    
    }

    	outfile << endl;

	std::array<std::vector<double>,5> v_energydep;
	std::array<std::vector<double>,5> v_x;
	std::array<std::vector<double>,5> v_y;

	int nhits_layers[5] = {0,0,0,0,0};

    for (int ev = 0;ev < nEvents; ev++) {
        tree->GetEntry(ev);


        int nHHits = Hhits->size();
       
        for (size_t j = 0; j < Hhits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = Hhits->at(j);
            double energy = hit->energyDeposit;
            double x = hit->position.x();
            double y = hit->position.y();
            double z = hit->position.z();
	

    	    int layer = DecodeLayer(z);	

	    v_x[layer].push_back(x);
	    v_y[layer].push_back(y);
	    v_energydep[layer].push_back(energy);
	    nhits_layers[layer]++;


	}

	outfile << ev << ",";

	for(int layerid = 0;layerid<5;layerid++){
		outfile << vec_av(v_x[layerid]) << "," << vec_av(v_y[layerid]) <<"," << vec_rms(v_x[layerid]) << "," << vec_rms(v_y[layerid]) << "," << nhits_layers[layerid] << vec_sum(v_energydep[layerid]) << "," << vec_rms(v_energydep[layerid]);
		
	}
	outfile << endl;
    }

outfile.close();
//    file->Close();
}
