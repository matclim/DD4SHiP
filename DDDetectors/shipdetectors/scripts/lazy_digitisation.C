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
#include <G4SystemOfUnits.hh> 

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


int DecodeLayerHCAL(double z){
	int layer = -9999;
	if(z > 500 && z < 1000) layer = 1;
	if(z > 1000 && z < 1200) layer = 2;
	if(z > 1200 && z < 1400) layer = 3;
	if(z > 1400 && z < 1600) layer = 4;
	if(z > 1600 && z < 2000) layer = 5;
    return layer;
}


int DecodeLayerSplitCALWide(double z){
	int layer = -9999;
	if(z > -600 && z < -500) layer = 2;
	if(z > -500 && z < -350) layer = 4;
	if(z > -350 && z < -200) layer = 6;
	if(z > -200 && z < 0) layer = 8;
	if(z > 0 && z < 200) layer = 10;
	if(z > 200 && z < 300) layer = 12;
	if(z > 300 && z < 450) layer = 14;
	if(z > 450 && z < 600) layer = 16;
    return layer;
}

int DecodeLayerSplitCALThin(double z){
	int layer = -9999;
	if(z > -600 && z < -500) layer = 1;
	if(z > -500 && z < -400) layer = 3;
	if(z > -400 && z < -300) layer = 5;
	if(z > -300 && z < 0) layer = 7;
	if(z > 0 && z < 150) layer = 9;
	if(z > 150 && z < 300) layer = 11;
	if(z > 300 && z < 400) layer = 13;
	if(z > 400 && z < 500) layer = 15;
	if(z > 500 && z < 600) layer = 17;
    return layer;
}

int DecodeBarWide(double x, double y, int layer){
	int bar = 0;
	
	if(layer%4 == 0){//Vertical bar
		if(x > -200 && x < -120) bar = 21;
		if(x > -120 && x < -50) bar = 22;
		if(x > -50 && x < 0) bar = 23;
		if(x > 0 && x < 60) bar = 24;
		if(x > 60 && x < 120) bar = 25;
		if(x > 120 && x < 200) bar = 26;
	}
	else{//Horizontal bar
		if(y > -200 && y < -120) bar = 11;
		if(y > -120 && y < -50) bar = 12;
		if(y > -50 && y < 0) bar = 13;
		if(y > 0 && y < 60) bar = 14;
		if(y > 60 && y < 120) bar = 15;
		if(y > 120 && y < 200) bar = 16;
	}
    return bar;
}


int DecodeBarThin(double x, int layer){
	int bar = 0;
	if(x > -100 && x < -90) bar = 1;
	if(x > -90 && x < -80) bar = 2;
	if(x > -80 && x < -70) bar = 3;
	if(x > -70 && x < -60) bar = 4;
	if(x > -60 && x < -50) bar = 5;
	if(x > -50 && x < -40) bar = 6;
	if(x > -40 && x < -30) bar = 7;
	if(x > -30 && x < -20) bar = 8;
	if(x > -20 && x < -10) bar = 9;
	if(x > -10 && x < 0) bar = 10;
	if(x > 0 && x < 10) bar = 11;
	if(x > 10 && x < 20) bar = 12;
	if(x > 20 && x < 30) bar = 13;
	if(x > 30 && x < 40) bar = 14;
	if(x > 40 && x < 50) bar = 15;
	if(x > 50 && x < 60) bar = 16;
	if(x > 60 && x < 70) bar = 17;
	if(x > 70 && x < 80) bar = 18;
	if(x > 80 && x < 90) bar = 19;
	if(x > 90 && x < 100) bar = 20;
    return bar;
}

int DecodeBarHCAL(double y, int layer){
	int bar = 0;
	if(y > -200 && y < -120) bar = 1;
	if(y > -120 && y < -70) bar = 2;
	if(y > -70 && y < 0) bar = 3;
	if(y > 0 && y < 70) bar = 4;
	if(y > 70 && y < 120) bar = 5;
	if(y > 120 && y < 200) bar = 6;
    return bar;
}

//Wide bar MPV: 1.71941 MeV
//Thin bar MPV: 1.69825 MeV

double calibration_wide = 1.71941;
double calibration_thin = 1.69825;

void lazy_digitisation(TString infile) {
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

    TFile* file = new TFile(infile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
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
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* whits = nullptr;
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* thits = nullptr;
    std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>* Hhits = nullptr;

    tree->SetBranchAddress("SplitCalWideBarHits", &whits);
    tree->SetBranchAddress("SplitCalThinBarHits", &thits);
    tree->SetBranchAddress("SHiPHCALHits", &Hhits);

    // ============================================================
    // 4. LOOP OVER EVENTS AND HITS
    // ============================================================
    int nEvents = tree->GetEntries();
    std::cout << "--- Starting analysis of " << nEvents << " events ---" << std::endl;
    std::cout << std::fixed << std::setprecision(5); // Set precision for coordinates/energy

    TH1F *h_wx = new TH1F("h_wx","h_wx",2000,-500,500);
    TH1F *h_wy = new TH1F("h_wy","h_wy",2000,-500,500);
    TH1F *h_wz = new TH1F("h_wz","h_wz",2000,-1000,1000);
   
    TH2F *h_wx_wz = new TH2F("h_wx_wz","h_wx_wz",100,-500,500,100,-1000,1000);
    TH2F *h_wy_wz= new TH2F("h_wy_wz","h_wy_wz",100,-500,500,100,-1000,1000);


    TH1F *h_tx = new TH1F("h_tx","h_tx",2000,-200,200);
    TH1F *h_ty = new TH1F("h_ty","h_ty",2000,-200,200);
    TH1F *h_tz = new TH1F("h_tz","h_tz",2000,-1000,1000);


    TH2F *h_tx_tz = new TH2F("h_tx_tz","h_tx_tz",500,-200,200,100,-1000,1000);
    TH2F *h_ty_tz= new TH2F("h_ty_tz","h_ty_tz",500,-200,200,100,-1000,1000);

    TH1F *h_hx = new TH1F("h_hx","h_hx",2000,-1000,1000);
    TH1F *h_hy = new TH1F("h_hy","h_hy",2000,-200,200);
    TH1F *h_hz = new TH1F("h_hz","h_hz",2000,0,2500);
   

    TH2F *h_hx_tz = new TH2F("h_hx_tz","h_hx_tz",500,-200,200,100,1000,2000);
    TH2F *h_hy_tz= new TH2F("h_hy_tz","h_hy_tz",500,-200,200,100,1000,2000);

    Long64_t event_id;
    std::vector<UShort_t> v_layers_wide;
    std::vector<UShort_t> v_bars_wide;
    std::vector<Double_t> v_n_mips_wide;

    std::vector<UShort_t> v_layers_thin;
    std::vector<UShort_t> v_bars_thin;
    std::vector<Double_t> v_n_mips_thin;
    
    std::vector<UShort_t> v_layers_hcal;
    std::vector<UShort_t> v_bars_hcal;
    std::vector<Double_t> v_n_mips_hcal;

    TTree *t_out = new TTree("lazy_digitised_events","Lazy diitisation of prototype calorimeter simulation events");

	
    t_out->Branch("event_id",&event_id);

    t_out->Branch("layers_wide",&v_layers_wide);
    t_out->Branch("bars_wide",&v_bars_wide);
    t_out->Branch("n_mips_wide",&v_n_mips_wide);

    t_out->Branch("layers_thin",&v_layers_thin);
    t_out->Branch("bars_thin",&v_bars_thin);
    t_out->Branch("n_mips_thin",&v_n_mips_thin);
    
    t_out->Branch("layers_hcal",&v_layers_hcal);
    t_out->Branch("bars_hcal",&v_bars_hcal);
    t_out->Branch("n_mips_hcal",&v_n_mips_hcal);


    TString ftitle = infile;

//    f_out = new TFile(ftitle.Remove(ftitle.Length()-5)+"_out.root");

    for (int ev = 0;ev < nEvents; ev++) {
        tree->GetEntry(ev);
	
	event_id = static_cast<Long64_t>(ev);

        int nwHits = whits->size();
        for (size_t j = 0; j < whits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = whits->at(j);
            double energy = hit->energyDeposit;
            double x = (hit->position.X());
            double y = hit->position.Y();
            double z = hit->position.Z();
	

	    h_wx->Fill(x);
	    h_wy->Fill(y);
	    h_wz->Fill(z);
	    
	    h_wx_wz->Fill(x,z);
	    h_wy_wz->Fill(y,z);
    	    UShort_t layer = DecodeLayerSplitCALWide(z);
    	    UShort_t bar   = DecodeBarWide(x,y,layer);
	    Double_t n_mips = energy / calibration_wide;

    	    v_layers_wide.push_back(layer);	    
    	    v_bars_wide.push_back(bar);	    
    	    v_n_mips_wide.push_back(n_mips);	    

	}
        
	int ntHits = thits->size();
        for (size_t j = 0; j < thits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = thits->at(j);
            double energy = hit->energyDeposit;
            double x = hit->position.X();
            double y = hit->position.Y();
            double z = hit->position.Z();

    	   cout << "Tz " << z << endl;	    

    	    //int layer = DecodeLayer(z);	

	    h_tx->Fill(x);
	    h_ty->Fill(y);
	    h_tz->Fill(z);
	    
	    h_tx_tz->Fill(x,z);
	    h_ty_tz->Fill(y,z);
    	    
	    UShort_t layer = DecodeLayerSplitCALThin(z);
    	    UShort_t bar   = DecodeBarThin(x,layer);
	    Double_t n_mips = energy / calibration_thin;

    	    v_layers_thin.push_back(layer);	    
    	    v_bars_thin.push_back(bar);	    
    	    v_n_mips_thin.push_back(n_mips);	    


	}
		

        int nHHits = Hhits->size();
        for (size_t j = 0; j < Hhits->size(); ++j) {
            dd4hep::sim::Geant4Calorimeter::Hit* hit = Hhits->at(j);
            double energy = hit->energyDeposit;
            double x = hit->position.X();
            double y = hit->position.Y();
            double z = hit->position.Z();
    	   
	    cout << "Hz " << z << endl;	    
	

    	    //int layer = DecodeLayer(z);	
	    h_hx->Fill(x);
	    h_hy->Fill(y);
	    h_hz->Fill(z);

	    h_hx_tz->Fill(x,z);
	    h_hy_tz->Fill(y,z);
	    
	    UShort_t layer = DecodeLayerSplitCALThin(z);
    	    UShort_t bar   = DecodeBarThin(y,layer);
	    Double_t n_mips = energy / calibration_wide;

    	    v_layers_hcal.push_back(layer);	    
    	    v_bars_hcal.push_back(bar);	    
    	    v_n_mips_hcal.push_back(n_mips);	    

	}

	t_out->Fill();

    	    v_layers_wide.clear();
    	    v_bars_wide.clear();
    	    v_n_mips_wide.clear();
    	    
	    v_layers_thin.clear();
    	    v_bars_thin.clear();
    	    v_n_mips_thin.clear();
    
    	    v_layers_hcal.clear();
    	    v_bars_hcal.clear();
    	    v_n_mips_hcal.clear();
    
    }


    //TCanvas *cwx = new TCanvas("cwx","cwx",800,600);
    //h_wx->Draw();
    //TCanvas *cwy = new TCanvas("cwy","cwy",800,600);
    //h_wy->Draw();
    //TCanvas *cwz = new TCanvas("cwz","cwz",800,600);
    //h_wz->Draw();

    //TCanvas *cwx_wz = new TCanvas("cwx_wz","cwx_wz",800,600);
    //h_wx_wz->Draw("COLZ");
    //TCanvas *cwy_wz = new TCanvas("cwy_wz","cwy_wz",800,600);
    //h_wy_wz->Draw("COLZ");

    //TCanvas *ctx = new TCanvas("ctx","ctx",800,600);
    //h_tx->Draw();
    //TCanvas *cty = new TCanvas("cty","cty",800,600);
    //h_ty->Draw();
    //TCanvas *ctz = new TCanvas("ctz","ctz",800,600);
    //h_tz->Draw();

    //TCanvas *ctx_tz = new TCanvas("ctx_tz","ctx_tz",800,600);
    //h_tx_tz->Draw("COLZ");
    //TCanvas *cty_tz = new TCanvas("cty_tz","cty_tz",800,600);
    //h_ty_tz->Draw("COLZ");


    //TCanvas *chx = new TCanvas("chx","chx",800,600);
    //h_hx->Draw();
    //TCanvas *chy = new TCanvas("chy","chy",800,600);
    //h_hy->Draw();
    //TCanvas *chz = new TCanvas("chz","chz",800,600);
    //h_hz->Draw();
    //
    //TCanvas *chx_tz = new TCanvas("chx_tz","chx_tz",800,600);
    //h_hx_tz->Draw("COLZ");
    //TCanvas *chy_tz = new TCanvas("chy_tz","chy_tz",800,600);
    //h_hy_tz->Draw("COLZ");


    ROOT::RDataFrame df(*t_out);


    df.Snapshot("df_lazy_sim_proto_DD4SHiP",ftitle.Remove(ftitle.Length()-5)+"_out.root");
}
