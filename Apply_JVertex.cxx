#include "JUNO_PMTs.h"

#include <numeric>
#include <TSpectrum.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include <vector>
#include <tuple>
#include <TTree.h>
#include <TTimeStamp.h>
#include "JVertex.h"

void execute(string filename, string treename, string outfilename, string PdfPath, float nEff) {
    // 1. Setup JVertex and PMTs
    string PMTFile = "/cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/data/Detector/Identifier/pmt_CDLPMT_latest.csv";
    
    // Instantiate vertex reconstructor once outside the loop
    JVertex vertex_reco(nEff, PdfPath, PMTFile);
    
    JUNO_PMTs PMTs_Pos; 
    PMTs_Pos.SetCdPmts(PMTFile);

    // 2. Setup Files and Trees
    TFile *fin = TFile::Open(filename.c_str(), "READ");
    TTree *intree = (TTree*)fin->Get(treename.c_str());
    TFile *file_out = TFile::Open(outfilename.c_str(), "RECREATE");
    TTree *out_tree = intree->CloneTree(0); 

    // 3. New Branches for Reconstruction Results
    float PromptX = -99999, PromptY = -99999, PromptZ = -99999;
    float DelayedX = -99999, DelayedY = -99999, DelayedZ = -99999;
    float PromptPE = 0, DelayedPE = 0;

    out_tree->Branch("PromptX", &PromptX, "PromptX/F");
    out_tree->Branch("PromptY", &PromptY, "PromptY/F");
    out_tree->Branch("PromptZ", &PromptZ, "PromptZ/F");
    out_tree->Branch("DelayedX", &DelayedX, "DelayedX/F");
    out_tree->Branch("DelayedY", &DelayedY, "DelayedY/F");
    out_tree->Branch("DelayedZ", &DelayedZ, "DelayedZ/F");
    out_tree->Branch("PromptPE", &PromptPE, "PromptPE/F");
    out_tree->Branch("DelayedPE", &DelayedPE, "DelayedPE/F");

    // 4. Connect input branches
    std::vector<float> *Time = nullptr, *Charge = nullptr, *CorrTime = nullptr;
    std::vector<int> *PMTID = nullptr, *PeakPositions = nullptr;
    float OECRecoX, OECRecoY, OECRecoZ;
    intree->SetBranchAddress("Time", &Time);
    intree->SetBranchAddress("Charge", &Charge);
    intree->SetBranchAddress("PMTID", &PMTID);
    intree->SetBranchAddress("RecoX", &OECRecoX);
    intree->SetBranchAddress("RecoY", &OECRecoY);
    intree->SetBranchAddress("RecoZ", &OECRecoZ);
	intree->SetBranchAddress("CorrTime", &CorrTime);
	intree->SetBranchAddress("PeakPositions", &PeakPositions);

    for (int i = 0; i < intree->GetEntries(); i++) {
        intree->GetEntry(i);
        
        // --- [Standard CorrTime Calculation & Peak Selection - Same as before] ---
        // (Assuming CorrTime_out and PeakPositions_out are calculated here)
        
        // Reset values for this event
        PromptX = PromptY = PromptZ = DelayedX = DelayedY = DelayedZ = -99999;
        PromptPE = DelayedPE = 0;

        if (PeakPositions->size() == 2) {
            // Temporary containers for filtered hits
            std::vector<float> p_Times, p_Charge, d_Times, d_Charge;
            std::vector<int> p_PMTID, d_PMTID;

            for (size_t j = 0; j < Charge->size(); j++) {
                float t = CorrTime ->at(j);
                
                // Filter Prompt Hits
                if (t < PeakPositions->at(0)*6 + 120 && t > PeakPositions->at(0)*6 - 60) {
                    p_Times.push_back(Time->at(j));
                    p_PMTID.push_back(PMTID->at(j));
                    p_Charge.push_back(Charge->at(j));
                    PromptPE += Charge->at(j);
                }
                // Filter Delayed Hits
                if (t < PeakPositions->at(1)*6 + 120 && t > PeakPositions->at(1)*6 - 60) {
                    d_Times.push_back(Time->at(j));
                    d_PMTID.push_back(PMTID->at(j));
                    d_Charge.push_back(Charge->at(j));
                    DelayedPE += Charge->at(j);
                }
            }

            // --- RECONSTRUCTION ---
            
            // 1. Reconstruct Prompt
            if (!p_Times.empty()) {
                vertex_reco.ChangeEvent(p_Times, p_PMTID, p_Charge);
                auto res = vertex_reco.GetEventPosition();
                PromptX = std::get<0>(res); PromptY = std::get<1>(res); PromptZ = std::get<2>(res);
            }

            // 2. Reconstruct Delayed
            if (!d_Times.empty()) {
                vertex_reco.ChangeEvent(d_Times, d_PMTID, d_Charge);
                auto res = vertex_reco.GetEventPosition();
                DelayedX = std::get<0>(res); DelayedY = std::get<1>(res); DelayedZ = std::get<2>(res);
            }
        }

        out_tree->Fill();
    }

    out_tree->Write();
    file_out->Close();
    fin->Close();
}


int main(int argc, char** argv) {

    string macro = argv[0];

    if(argc<5 || argc > 6) {
            cout << "\n     USAGE:  "<< macro << " Input_Rootfile  In_Tree_Name Out_File_Name JVertex_PDF_Path [nEff_for_JVertex] \n" << endl;
            return 1;
    }

    string filename = argv[1];
	string intreename = argv[2];
	string outfilename = argv[3];
	string PdfPath = argv[4];

	float nEff;

	if (argc == 6) {
		nEff = stof(argv[5]);
	} else {
		nEff = 1.60;
	}

	execute (filename, intreename, outfilename, PdfPath, nEff);

    return 0;
}

