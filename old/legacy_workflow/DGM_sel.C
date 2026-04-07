#include <sys/stat.h>
#include <iostream>
#include <filesystem>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooFit.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm> // For std::lower_bound
using namespace std;
using namespace RooFit;


namespace fs = std::filesystem;
/////////////////////////////////////
////// User defined variables ////// 
// Specify data type
std::string datatype = "DATA";
// Specify Muon type
std::string Muon_type = "DSA";
////////////////////////////////////


// Global constants
// pT
const std::array<double, 7> PT_BINS = { 20, 40, 55, 70, 100, 150, 1000 };
float binSize[6] = { 0.0 };
float binCenter[6] = { 0.0 };
// |dZ|
const std::array<double, 7> DZ_BINS = { 1, 10, 20, 30, 45, 60, 150 };
float DZ_binSize[6] = { 0.0 };
float DZ_binCenter[6] = { 0.0 };
// |dXY|
const std::array<double, 6> DXY_BINS = { 1, 10, 20, 30, 40, 80 };
float DXY_binSize[5] = { 0.0 };
float DXY_binCenter[5] = { 0.0 };


// Function declarations
void InitializeHistograms(std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F*& histPtTag, TH1F*& histEtaTag, TH1F*& histPhiTag, TH1F*& histChargeTag, TH1F*& histPtProbe, TH1F*& histEtaProbe, TH1F*& histPhiProbe, TH1F*& histChargeProbe, TH1F*& histPtDGL, const std::string& Muon);

void ProcessFiles(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag, TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL);

void ProcessFiles_DSA(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag, TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDSA);

void FitAndDrawHistograms(const std::vector<TH1F*>& ptHistograms, TFile* outputFile, const std::vector<double>& _BINS, const std::string& bin_type, float* binCenter, float* binSize, const char* X_tile, const std::string& Muon, const std::string& data_type);

void DrawControlPlots(TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag, TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL, TFile* outputFile, const std::string& Muon, const std::string& data_type);

// Main Function
void DGM_sel() {


    // Create a ROOT file to save histograms
    TFile* outputFile = new TFile(("Cosmics_muons_" + datatype + "_" + Muon_type + ".root").c_str(), "RECREATE");

    // List of directories containing ROOT files
    std::vector<std::string> directories;
    if (datatype == "DATA") {
        directories = {
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023C_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171725/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023D_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171752/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023E_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171805/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023F_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171822/0000"
        };
    }
    else if (datatype == "MC") {
        directories = {
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/UndergroundCosmiHPLooseMu_bottomPhiFilter/CosmicsAnalysis_Run2022_MC_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_MC_2025_customTTrig/250217_181515/0000"
        };
    }
    else {
        std::cerr << "Invalid datatype: " << datatype << std::endl;
        return;
    }

    // Vector to store all ROOT files from all directories
    std::vector<std::string> files;

    // Loop over directories to collect ROOT files
    for (const auto& dir : directories) {
        if (fs::is_directory(dir)) {
            for (const auto& entry : fs::directory_iterator(dir)) {
                if (entry.path().extension() == ".root") {
                    files.push_back(entry.path().string());
                    std::cout << "Found file: " << entry.path() << std::endl;
                }
            }
        }
        else {
            std::cerr << "Invalid directory: " << dir << std::endl;
        }
    }

    // Histograms for pT bins
    std::vector<TH1F*> ptHistograms;
    std::vector<TH1F*> dzHistograms;
    std::vector<TH1F*> dxyHistograms;
    TH1F* histPtTag, * histEtaTag, * histPhiTag, * histChargeTag;
    TH1F* histPtProbe, * histEtaProbe, * histPhiProbe, * histChargeProbe;
    TH1F* histPtDGL;
    TH1F* Tot_pthist;

    for (int n = 0; n < PT_BINS.size() - 1; n++) {
        binSize[n] = (PT_BINS[n + 1] - PT_BINS[n]) / 2;
        binCenter[n] = binSize[n] + PT_BINS[n];

    }
    for (int n = 0; n < DZ_BINS.size() - 1; n++) {
        DZ_binSize[n] = (DZ_BINS[n + 1] - DZ_BINS[n]) / 2;
        DZ_binCenter[n] = DZ_binSize[n] + DZ_BINS[n];

    }

    for (int n = 0; n < DXY_BINS.size() - 1; n++) {
        DXY_binSize[n] = (DXY_BINS[n + 1] - DXY_BINS[n]) / 2;
        DXY_binCenter[n] = DXY_binSize[n] + DXY_BINS[n];

    }


    // Initialize histograms
    std::vector<double> pt_bins(std::begin(PT_BINS), std::end(PT_BINS));
    std::vector<double> dz_bins(std::begin(DZ_BINS), std::end(DZ_BINS));
    std::vector<double> dxy_bins(std::begin(DXY_BINS), std::end(DXY_BINS));

    InitializeHistograms(ptHistograms, dzHistograms, dxyHistograms, Tot_pthist, histPtTag, histEtaTag, histPhiTag, histChargeTag, histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL, Muon_type);

    // Process files
    if (Muon_type == "DGL") {
        ProcessFiles(files, ptHistograms, dzHistograms, dxyHistograms, Tot_pthist, histPtTag, histEtaTag, histPhiTag, histChargeTag, histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL);
    }
    else if (Muon_type == "DSA") {
        ProcessFiles_DSA(files, ptHistograms, dzHistograms, dxyHistograms, Tot_pthist, histPtTag, histEtaTag, histPhiTag, histChargeTag, histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL);

    }

    // Fit and draw histograms
    FitAndDrawHistograms(ptHistograms, outputFile, pt_bins, "pt", binCenter, binSize, "p_{T} #mu_{ref} [GeV] ", Muon_type, datatype);
    //FitAndDrawHistograms(dzHistograms, outputFile, dz_bins, "dz", DZ_binCenter, DZ_binSize, "|dz|", Muon_type, datatype);
    //FitAndDrawHistograms(dxyHistograms, outputFile, dxy_bins, "dxy", DXY_binCenter, DXY_binSize, "|dxy|", Muon_type, datatype);

    // Draw control plots
    DrawControlPlots(histPtTag, histEtaTag, histPhiTag, histChargeTag, histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL, outputFile, Muon_type, datatype);

    // Clean up
    for (auto hist : ptHistograms) delete hist;
    for (auto hist : dzHistograms) delete hist;
    for (auto hist : dxyHistograms) delete hist;
    delete histPtTag;
    delete histPtProbe;

    // Close the output file
    outputFile->Close();
    delete outputFile;

    std::cout << "Finished processing and histogram saving." << std::endl;
}

void InitializeHistograms(std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F*& histPtTag, TH1F*& histEtaTag, TH1F*& histPhiTag, TH1F*& histChargeTag,
    TH1F*& histPtProbe, TH1F*& histEtaProbe, TH1F*& histPhiProbe, TH1F*& histChargeProbe, TH1F*& histPtDGL, const std::string& Muon) {
    float min_bin = 0.0;
    float max_bin = 0.0;
    if (Muon == "DGL") {
        min_bin = -0.3;
        max_bin = 0.3;

    }
    else if (Muon == "DSA") {
        min_bin = -5.0;
        max_bin = 5.0;

    }

    // Initialize histograms for pT each bin
    for (size_t i = 0; i < PT_BINS.size() - 1; ++i) {
        ptHistograms.push_back(new TH1F(Form("pt_%g-%g", PT_BINS[i], PT_BINS[i + 1]),
            Form("pT %g-%g", PT_BINS[i], PT_BINS[i + 1]), 100, min_bin, max_bin));
    }

    // Initialize histograms for dz each bin
    for (size_t i = 0; i < DZ_BINS.size() - 1; ++i) {
        dzHistograms.push_back(new TH1F(Form("|dZ|_%g-%g", DZ_BINS[i], DZ_BINS[i + 1]),
            Form("|dz| %g-%g", DZ_BINS[i], DZ_BINS[i + 1]), 100, min_bin, max_bin));
    }

    // Initialize histograms for dz each bin
    for (size_t i = 0; i < DXY_BINS.size() - 1; ++i) {
        dxyHistograms.push_back(new TH1F(Form("|dXY|_%g-%g", DXY_BINS[i], DXY_BINS[i + 1]),
            Form("|dXY| %g-%g", DXY_BINS[i], DXY_BINS[i + 1]), 100, min_bin, max_bin));
    }

    // Initialize other histograms
    Tot_pthist = new TH1F("Tot_pthist", "q/pT Residual distribution", 80, -0.3, 0.3);
    histPtTag = new TH1F("hist_pt_tag", "Tag muon pT", 350, 0, 300);
    histEtaTag = new TH1F("hist_eta_tag", "Tag muon #eta", 60, -3, 3);
    histPhiTag = new TH1F("hist_phi_tag", "Tag muon #phi", 60, -3, -3);
    histChargeTag = new TH1F("hist_charge_tag", "Tag muon q", 10, -4, 4);

    histPtProbe = new TH1F("hist_pt_probe", "Probe muon pT", 350, 0, 300);
    histEtaProbe = new TH1F("hist_eta_probe", "Probe muon #eta", 60, -3, 3);
    histPhiProbe = new TH1F("hist_phi_probe", "Probe muon #phi", 60, -3, -3);
    histChargeProbe = new TH1F("hist_charge_probe", "Probe muon q", 10, -4, 4);

    histPtDGL = new TH1F(("hist_pt_" + Muon).c_str(), (Muon + " muon pT").c_str(), 350, 0, 300);
    std::cout << "Histograms initialized" << std::endl;
}

void ProcessFiles(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
    TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL) {
    // Loop through each file
    for (const auto& file : files) {
        std::cout << "Reading file: " << file << std::endl;
        TFile* f = TFile::Open(file.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << file << std::endl;
            continue;
        }

        // Get the TTree
        TTree* tree = dynamic_cast<TTree*>(f->Get("Events"));
        if (!tree) {
            std::cerr << "Error: TTree 'Events' not found!" << std::endl;
            f->Close();
            delete f;
            continue;
        }

        // Variables to hold branch data
        const int maxMuons = 200; // Adjust this to the maximum expected size of ndmu
        int ndmu;
        int event;
        int dmuNumberOfChambersCSCorDT[maxMuons];
        float dmuDglPt[maxMuons];
        float dmuDglEta[maxMuons];
        float dmuDglPhi[maxMuons];
        float dmuDglDz[maxMuons];
        float dmuDglDxy[maxMuons];
        bool dmuDglPassTagID[maxMuons];
        bool dmuDglHasProbe[maxMuons];
        int dmuDglProbeID[maxMuons];
        int dmuIsDGL[maxMuons];
        float dmuDglCharge[maxMuons];
        bool hltL2Mu10NoVertexNoBPTX3BX = false;

        // Link branches
        tree->SetBranchAddress("ndmu", &ndmu);
        tree->SetBranchAddress("event", &event);
        tree->SetBranchAddress("dmu_numberOfChambersCSCorDT", dmuNumberOfChambersCSCorDT);
        tree->SetBranchAddress("dmu_dgl_pt", dmuDglPt);
        tree->SetBranchAddress("dmu_dgl_eta", dmuDglEta);
        tree->SetBranchAddress("dmu_dgl_phi", dmuDglPhi);
        tree->SetBranchAddress("dmu_dgl_dz", dmuDglDz);
        tree->SetBranchAddress("dmu_dgl_dxy", dmuDglDxy);
        tree->SetBranchAddress("dmu_dgl_charge", dmuDglCharge);
        tree->SetBranchAddress("dmu_dgl_passTagID", dmuDglPassTagID);
        tree->SetBranchAddress("dmu_dgl_probeID", dmuDglProbeID);
        tree->SetBranchAddress("dmu_dgl_hasProbe", dmuDglHasProbe);
        tree->SetBranchAddress("dmu_isDGL", dmuIsDGL);
        tree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &hltL2Mu10NoVertexNoBPTX3BX);

        // Loop over all entries
        int nEntries = tree->GetEntries();
        for (int entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);

            // Only select events that pass the HLT path
            if (datatype == "MC") hltL2Mu10NoVertexNoBPTX3BX = true;
            if (!hltL2Mu10NoVertexNoBPTX3BX) continue;

            // Tag variables
            float tagPt = 0.0;
            float tagEta = 0.0;
            float tagPhi = 0.0;
            float tagDz = 0.0;
            float tagDxy = 0.0;
            float tagCharge = 0.0;

            // Probe variables
            float probePt = 0.0;
            float probeEta = 0.0;
            float probePhi = 0.0;
            float probeDz = 0.0;
            float probeDxy = 0.0;
            float probeCharge = 0.0;

            // Select events with 2 muons or more
            if (ndmu < 2) continue;

            // Loop over the number of muons (ndmu) and select the Tag muon
            for (int i = 0; i < ndmu; i++) {
                histPtDGL->Fill(dmuDglPt[i]);

                if (!dmuDglPassTagID[i]) continue;

                if (dmuDglPassTagID[i] && dmuIsDGL[i] == 1) {
                    // Fill tag histograms
                    tagPt = dmuDglPt[i];
                    tagEta = dmuDglEta[i];
                    tagPhi = dmuDglPhi[i];
                    tagDz = abs(dmuDglDz[i]);
                    tagDxy = abs(dmuDglDxy[i]);
                    tagCharge = dmuDglCharge[i];
                    histPtTag->Fill(dmuDglPt[i]);
                    histEtaTag->Fill(dmuDglEta[i]);
                    histPhiTag->Fill(dmuDglPhi[i]);
                    histChargeTag->Fill(dmuDglCharge[i]);

                    // Probe selection
                    if (dmuDglHasProbe[i]) {
                        int probeIndex = dmuDglProbeID[i];
                        probePt = dmuDglPt[probeIndex];
                        probeEta = dmuDglEta[probeIndex];
                        probePhi = dmuDglPhi[probeIndex];
                        probeDz = dmuDglDz[probeIndex];
                        probeDxy = dmuDglDxy[probeIndex];
                        probeCharge = dmuDglCharge[probeIndex];
                        histPtProbe->Fill(dmuDglPt[probeIndex]);
                        histEtaProbe->Fill(dmuDglEta[probeIndex]);
                        histPhiProbe->Fill(dmuDglPhi[probeIndex]);
                        histChargeProbe->Fill(dmuDglCharge[probeIndex]);

                        // Calculate the resolution
                        float invUppt = std::abs(dmuDglCharge[probeIndex] / dmuDglPt[probeIndex]);
                        float invDownpt = std::abs(dmuDglCharge[i] / dmuDglPt[i]);
                        double resolution = (invUppt - invDownpt) / (std::sqrt(2) * invDownpt);

                        //fill q/pT residual hist
                        Tot_pthist->Fill(resolution);

                        // Fill pT histograms
                        auto it = std::lower_bound(PT_BINS.begin(), PT_BINS.end(), tagPt);
                        size_t index = std::distance(PT_BINS.begin(), it);

                        if (index > 0 && index < PT_BINS.size()) {
                            ptHistograms[index - 1]->Fill(resolution);
                        }
                        /*for (size_t j = 0; j < PT_BINS.size() - 1; ++j) {
                            if (PT_BINS[j] < tagPt && tagPt < PT_BINS[j + 1]) {
                                ptHistograms[j]->Fill(resolution);
                                break;  // Exit loop once the correct bin is found
                            }
                        }*/

                        // Fill dz histograms
                        auto it_dz = std::lower_bound(DZ_BINS.begin(), DZ_BINS.end(), tagDz);
                        size_t index_dz = std::distance(DZ_BINS.begin(), it_dz);

                        if (index_dz > 0 && index_dz < DZ_BINS.size()) {
                            dzHistograms[index_dz - 1]->Fill(resolution);
                        }
                        /*for (size_t j = 0; j < DZ_BINS.size() - 1; ++j) {
                            if (DZ_BINS[j] < tagDz && tagDz < DZ_BINS[j + 1]) {
                                dzHistograms[j]->Fill(resolution);
                                break;  // Exit loop once the correct bin is found
                            }
                        }*/

                        // Fill dxy histograms
                        auto it_dyx = std::lower_bound(DXY_BINS.begin(), DXY_BINS.end(), tagDxy);
                        size_t index_dyx = std::distance(DXY_BINS.begin(), it_dyx);

                        if (index_dyx > 0 && index_dyx < DXY_BINS.size()) {
                            dxyHistograms[index_dyx - 1]->Fill(resolution);
                        }
                        /*for (size_t j = 0; j < DXY_BINS.size() - 1; ++j) {
                            if (DXY_BINS[j] < tagDxy && tagDxy < DXY_BINS[j + 1]) {
                                dxyHistograms[j]->Fill(resolution);
                                break;  // Exit loop once the correct bin is found
                            }
                        }*/
                    }
                }
            }
        }

        // Close the file
        f->Close();
        delete f;
    }
    std::cout << "Reading files " << std::endl;
    //ptHistograms.push_back(Tot_pthist);
}


void ProcessFiles_DSA(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, std::vector<TH1F*>& dzHistograms, std::vector<TH1F*>& dxyHistograms, TH1F*& Tot_pthist, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
    TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDSA) {
    // Loop through each file
    for (const auto& file : files) {
        std::cout << "Reading file: " << file << std::endl;
        TFile* f = TFile::Open(file.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << file << std::endl;
            continue;
        }

        // Get the TTree
        TTree* tree = dynamic_cast<TTree*>(f->Get("Events"));
        if (!tree) {
            std::cerr << "Error: TTree 'Events' not found!" << std::endl;
            f->Close();
            delete f;
            continue;
        }

        // Variables to hold branch data
        const int maxMuons = 200; // Adjust this to the maximum expected size of ndmu
        int ndmu;
        int event;
        int dmuNumberOfChambersCSCorDT[maxMuons];
        float dmuDsaPt[maxMuons];
        float dmuDsaEta[maxMuons];
        float dmuDsaPhi[maxMuons];
        float dmuDsaDz[maxMuons];
        float dmuDsaDxy[maxMuons];
        bool dmuDsaPassTagID[maxMuons];
        bool dmuDsaHasProbe[maxMuons];
        int dmuDsaProbeID[maxMuons];
        int dmuIsDSA[maxMuons];
        float dmuDsaCharge[maxMuons];
        bool hltL2Mu10NoVertexNoBPTX3BX = false;

        // Link branches
        tree->SetBranchAddress("ndmu", &ndmu);
        tree->SetBranchAddress("event", &event);
        tree->SetBranchAddress("dmu_numberOfChambersCSCorDT", dmuNumberOfChambersCSCorDT);
        tree->SetBranchAddress("dmu_dsa_pt", dmuDsaPt);
        tree->SetBranchAddress("dmu_dsa_eta", dmuDsaEta);
        tree->SetBranchAddress("dmu_dsa_phi", dmuDsaPhi);
        tree->SetBranchAddress("dmu_dsa_dz", dmuDsaDz);
        tree->SetBranchAddress("dmu_dsa_dxy", dmuDsaDxy);
        tree->SetBranchAddress("dmu_dsa_charge", dmuDsaCharge);
        tree->SetBranchAddress("dmu_dsa_passTagID", dmuDsaPassTagID);
        tree->SetBranchAddress("dmu_dsa_probeID", dmuDsaProbeID);
        tree->SetBranchAddress("dmu_dsa_hasProbe", dmuDsaHasProbe);
        tree->SetBranchAddress("dmu_isDSA", dmuIsDSA);
        tree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &hltL2Mu10NoVertexNoBPTX3BX);

        // Loop over all entries
        int nEntries = tree->GetEntries();
        for (int entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);

            // Only select events that pass the HLT path
            if (datatype == "MC") hltL2Mu10NoVertexNoBPTX3BX = true;
            if (!hltL2Mu10NoVertexNoBPTX3BX) continue;

            // Tag variables
            float tagPt = 0.0;
            float tagEta = 0.0;
            float tagPhi = 0.0;
            float tagDz = 0.0;
            float tagDxy = 0.0;
            float tagCharge = 0.0;

            // Probe variables
            float probePt = 0.0;
            float probeEta = 0.0;
            float probePhi = 0.0;
            float probeDz = 0.0;
            float probeDxy = 0.0;
            float probeCharge = 0.0;

            // Select events with 2 muons or more
            if (ndmu < 2) continue;

            // Loop over the number of muons (ndmu) and select the Tag muon
            for (int i = 0; i < ndmu; i++) {
                histPtDSA->Fill(dmuDsaPt[i]);

                if (dmuDsaPt[i] < 12.5) continue;

                if (!dmuDsaPassTagID[i]) continue;

                if (dmuDsaPassTagID[i] && dmuIsDSA[i] == 1) {
                    // Fill tag histograms
                    tagPt = dmuDsaPt[i];
                    tagEta = dmuDsaEta[i];
                    tagPhi = dmuDsaPhi[i];
                    tagDz = abs(dmuDsaDz[i]);
                    tagDxy = abs(dmuDsaDxy[i]);
                    tagCharge = dmuDsaCharge[i];
                    histPtTag->Fill(dmuDsaPt[i]);
                    histEtaTag->Fill(dmuDsaEta[i]);
                    histPhiTag->Fill(dmuDsaPhi[i]);
                    histChargeTag->Fill(dmuDsaCharge[i]);

                    // Probe selection
                    if (dmuDsaHasProbe[i]) {
                        int probeIndex = dmuDsaProbeID[i];
                        if (dmuDsaPt[probeIndex] < 12.5) continue;
                        probePt = dmuDsaPt[probeIndex];
                        probeEta = dmuDsaEta[probeIndex];
                        probePhi = dmuDsaPhi[probeIndex];
                        probeDz = dmuDsaDz[probeIndex];
                        probeDxy = dmuDsaDxy[probeIndex];
                        probeCharge = dmuDsaCharge[probeIndex];
                        histPtProbe->Fill(dmuDsaPt[probeIndex]);
                        histEtaProbe->Fill(dmuDsaEta[probeIndex]);
                        histPhiProbe->Fill(dmuDsaPhi[probeIndex]);
                        histChargeProbe->Fill(dmuDsaCharge[probeIndex]);



                        // Calculate the resolution
                        float invUppt = std::abs(dmuDsaCharge[probeIndex] / dmuDsaPt[probeIndex]);
                        float invDownpt = std::abs(dmuDsaCharge[i] / dmuDsaPt[i]);
                        double resolution = (invUppt - invDownpt) / (std::sqrt(2) * invDownpt);

                        //fill q/pT residual hist
                        Tot_pthist->Fill(resolution);

                        // Fill pT histograms
                        auto it = std::lower_bound(PT_BINS.begin(), PT_BINS.end(), tagPt);
                        size_t index = std::distance(PT_BINS.begin(), it);

                        if (index > 0 && index < PT_BINS.size()) {
                            ptHistograms[index - 1]->Fill(resolution);
                        }

                        // Fill dz histograms
                        auto it_dz = std::lower_bound(DZ_BINS.begin(), DZ_BINS.end(), tagDz);
                        size_t index_dz = std::distance(DZ_BINS.begin(), it_dz);

                        if (index_dz > 0 && index_dz < DZ_BINS.size()) {
                            dzHistograms[index_dz - 1]->Fill(resolution);
                        }

                        // Fill dxy histograms
                        auto it_dyx = std::lower_bound(DXY_BINS.begin(), DXY_BINS.end(), tagDxy);
                        size_t index_dyx = std::distance(DXY_BINS.begin(), it_dyx);

                        if (index_dyx > 0 && index_dyx < DXY_BINS.size()) {
                            dxyHistograms[index_dyx - 1]->Fill(resolution);

                        }
                    }
                }
            }
        }

        // Close the file
        f->Close();
        delete f;
    }
    std::cout << "Reading files " << std::endl;
    //ptHistograms.push_back(Tot_pthist);
}

void FitAndDrawHistograms(const std::vector<TH1F*>& ptHistograms, TFile* outputFile, const std::vector<double>& _BINS, const std::string& bin_type, float* binCenter, float* binSize, const char* X_tile, const std::string& Muon, const std::string& data_type) {
    // Mean and Sigma arrays
    int nd = ptHistograms.size();
    std::vector<float> mean(nd, 0.0);
    std::vector<float> meanErr(nd, 0.0);
    std::vector<float> sigma(nd, 0.0);
    std::vector<float> sigmaErr(nd, 0.0);

    outputFile->cd();
    std::cout << "Making plots " << std::endl;
    // Create a canvas for all plots together
    TCanvas* canvas_all = new TCanvas(("All_plots_" + bin_type).c_str(), (bin_type + " Histograms").c_str(), 1500, 1000);
    canvas_all->Divide(3, 2);  // Adjust layout as needed (4 columns, 2 rows)
    canvas_all->SetGrid();
    float min_bin = 0.0;
    float max_bin = 0.0;
    if (Muon == "DGL") {
        min_bin = -0.3;
        max_bin = 0.3;

    }
    else if (Muon == "DSA") {
        min_bin = -5.0;
        max_bin = 5.0;
    }
    for (size_t i = 0; i < ptHistograms.size(); ++i) {
        std::string canvas_name = "c_" + std::string(bin_type) + "_" + std::to_string(i);
        std::string histogram_title = std::string(bin_type) + " Histogram " + std::to_string(_BINS[i]) + "-" + std::to_string(_BINS[i + 1]);

        TCanvas* canvas = new TCanvas(canvas_name.c_str(), histogram_title.c_str(), 800, 600);
        canvas->SetGrid();

        // Fit the histogram with a Gaussian
        TF1* gaussFit = new TF1("gaussFit", "gaus", min_bin, max_bin);
        ptHistograms[i]->Fit(gaussFit, "SR");

        // Extract fit parameters
        double chi2 = gaussFit->GetChisquare();
        int ndf = gaussFit->GetNDF();
        double chi2PerNdf = chi2 / ndf;
        double meanVal = gaussFit->GetParameter(1);
        double stdDev = gaussFit->GetParameter(2);

        mean[i] = meanVal;
        meanErr[i] = gaussFit->GetParError(1);
        sigma[i] = stdDev;
        sigmaErr[i] = gaussFit->GetParError(2);

        // Print parameters
        std::cout << "Histogram: " << bin_type << _BINS[i] << "-" << _BINS[i + 1]
            << ", Mean: " << meanVal << ", Std Dev: " << stdDev << std::endl;

        // Draw histogram and fit
        ptHistograms[i]->Draw();
        ptHistograms[i]->GetYaxis()->SetTitle("Events");
        ptHistograms[i]->GetXaxis()->SetTitle("q/p_{T} Residuals");
        gaussFit->Draw("SAME");

        // Add legend with parameters
        TLegend* legend = new TLegend(0.6, 0.5, 0.9, 0.75);
        legend->SetTextSize(0.04);
        legend->SetFillStyle(0);
        legend->AddEntry((TObject*)nullptr, Form("Mean = %.4f", meanVal), "");
        legend->AddEntry((TObject*)nullptr, Form("#sigma = %.4f", stdDev), "");
        legend->AddEntry((TObject*)nullptr, Form("#chi^2 = %.4f", chi2), "");
        legend->AddEntry((TObject*)nullptr, Form("NDF = %d", ndf), "");
        legend->AddEntry((TObject*)nullptr, Form("#chi^2/ndf = %.4f", chi2PerNdf), "");
        legend->SetBorderSize(0);
        legend->Draw();

        ptHistograms[i]->Write();

        // Save canvas
        std::string save_path = Form("%s/%s/%s_Plots/dmu_%s_%s_%.0f-%.0f.png", data_type.c_str(), Muon.c_str(), bin_type.c_str(), bin_type.c_str(), Muon.c_str(), _BINS[i], _BINS[i + 1]);

        canvas->SaveAs(save_path.c_str());

        // Add the histogram to the combined canvas
        canvas_all->cd(i + 1);
        ptHistograms[i]->Draw();
        gaussFit->Draw("SAME");

        delete gaussFit;
        delete canvas;
    }

    // Draw all canvas with all plots 
    canvas_all->cd();
    canvas_all->SaveAs(Form("%s/%s/%s_Plots/dmu_%s_%s_All_fits.png", data_type.c_str(), Muon.c_str(), bin_type.c_str(),  Muon.c_str(), bin_type.c_str()));

    // Draw mean and sigma plots
    TCanvas* cc1 = new TCanvas("cc1", "cc1", 1000, 750);
    cc1->cd();
    cc1->SetGrid();
    cc1->SetLogx();
    cc1->SetLeftMargin(0.24);


    auto meanGraph = new TGraphErrors(nd, binCenter, mean.data(), binSize, meanErr.data());
    meanGraph->SetTitle("mean_total");
    meanGraph->SetLineWidth(2);
    meanGraph->SetMarkerStyle(8);
    meanGraph->GetXaxis()->SetTitle(X_tile);
    meanGraph->GetYaxis()->SetTitle("Mean of q/p_{T} relative residual ");
    meanGraph->Draw("A*");
    meanGraph->Write((bin_type + "mean_total").c_str());
    //cc1->SaveAs((bin_type + "mean_total_.png").c_str());
    delete cc1;
    TCanvas* cc2 = new TCanvas("cc2", "cc2", 1000, 750);
    cc2->cd();
    cc2->SetGrid();
    cc2->SetLogx();
    cc2->SetLeftMargin(0.24);

    auto sigmaGraph = new TGraphErrors(nd, binCenter, sigma.data(), binSize, sigmaErr.data());
    sigmaGraph->SetTitle("Sigma_total");
    sigmaGraph->SetLineWidth(2);
    sigmaGraph->SetMarkerStyle(8);
    sigmaGraph->GetXaxis()->SetTitle(X_tile);
    sigmaGraph->GetYaxis()->SetTitle("#sigma of q/p_{T} relative residual ");
    sigmaGraph->Draw("A*");
    sigmaGraph->Write((bin_type + "Sigma_total").c_str());
    //cc2->SaveAs((bin_type + "Sigma_total_.png").c_str());
    delete cc1;
}


void DrawControlPlots(TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
    TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL, TFile* outputFile, const std::string& Muon, const std::string& data_type) {
    // Draw control plots for Tag muons
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    histPtTag->Draw();
    histPtTag->GetYaxis()->SetTitle("Events");
    histPtTag->GetXaxis()->SetTitle("p_{T}");
    histPtTag->Write();
    c1->SaveAs((data_type + "/"+ Muon + "/Control_plots/hist_pt_tag.png").c_str());

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    histEtaTag->Draw();
    histEtaTag->GetYaxis()->SetTitle("Events");
    histEtaTag->GetXaxis()->SetTitle("#eta");
    histEtaTag->Write();
    c2->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_eta_tag.png").c_str());

    TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    histPhiTag->Draw();
    histPhiTag->GetYaxis()->SetTitle("Events");
    histPhiTag->GetXaxis()->SetTitle("#phi");
    histPhiTag->Write();
    c3->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_phi_tag.png").c_str());

    TCanvas* c63 = new TCanvas("c63", "c63", 800, 600);
    histChargeTag->Draw();
    histChargeTag->GetYaxis()->SetTitle("Events");
    histChargeTag->GetXaxis()->SetTitle("charge");
    histChargeTag->Write();
    c63->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_charge_tag.png").c_str());

    // Draw control plots for Probe muons
    TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
    histPtProbe->Draw();
    histPtProbe->GetYaxis()->SetTitle("Events");
    histPtProbe->GetXaxis()->SetTitle("p_{T}");
    histPtProbe->Write();
    c4->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_pt_probe.png").c_str());

    TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
    histEtaProbe->Draw();
    histEtaProbe->GetYaxis()->SetTitle("Events");
    histEtaProbe->GetXaxis()->SetTitle("#eta");
    histEtaProbe->Write();
    c5->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_eta_probe.png").c_str());

    TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
    histPhiProbe->Draw();
    histPhiProbe->GetYaxis()->SetTitle("Events");
    histPhiProbe->GetXaxis()->SetTitle("#phi");
    histPhiProbe->Write();
    c6->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_phi_probe.png").c_str());

    TCanvas* c66 = new TCanvas("c66", "c66", 800, 600);
    histChargeProbe->Draw();
    histChargeProbe->GetYaxis()->SetTitle("Events");
    histChargeProbe->GetXaxis()->SetTitle("charge");
    histChargeProbe->Write();
    c66->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_charge_probe.png").c_str());

    TCanvas* c7 = new TCanvas("c7", "c7", 800, 600);
    histPtDGL->Draw();
    histPtDGL->Write();
    c7->SaveAs((data_type + "/" + Muon + "/Control_plots/hist_pT_DGL.png").c_str());
}
