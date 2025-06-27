#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace Pythia8;
using namespace fastjet;

// Function to trace back a particle to its hard process ancestor
int traceToHardProcess(const Particle& particle, const Event& event) {
    int current = particle.mother1();
    while (current > 0) {
        const Particle& mother = event[current];
        // Check if this is a parton from the hard process
        if (mother.statusAbs() == 23 || mother.statusAbs() == 21) {
            return mother.id();
        }
        current = mother.mother1();
    }
    return 0; // Return 0 if no hard process ancestor found
}

// Get the most common hard process ancestor among jet constituents
int getJetPartonID(const PseudoJet& jet, const Event& event) {
    std::map<int, int> idCounts;
    
    // Loop over jet constituents
    std::vector<PseudoJet> constituents = jet.constituents();
    for (const PseudoJet& constituent : constituents) {
        // Get the Pythia particle index (assuming user_index was set)
        int idx = constituent.user_index();
        if (idx < 0 || idx >= event.size()) continue;
        
        const Particle& part = event[idx];
        int hardID = traceToHardProcess(part, event);
        if (hardID != 0) {
            idCounts[hardID]++;
        }
    }
    
    // Find the most common hard process ID
    if (idCounts.empty()) return 0;
    
    auto maxEntry = std::max_element(idCounts.begin(), idCounts.end(),
        [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.second < b.second;
        });
    
    return maxEntry->first;
}

int main(int argc, char* argv[]) {
    // Parse random seed and lambda value from command line
    int seed = 0;
    int lambdaT = 10000;
    if (argc > 1) seed = std::stoi(argv[1]);
    if (argc > 2) lambdaT = std::stoi(argv[2]);

    // Initialize Pythia
    Pythia pythia;
    // pythia.readString("HardQCD:all = on"); HardQCD is already included in led framework
    pythia.readString("ExtraDimensionsLED:dijets = on");

    pythia.readString("ExtraDimensionsLED:n = 6");
    pythia.readString("ExtraDimensionsLED:MD = 2000"); //Default value - fundamental scale of gravity in D = 4 + n dimensions
    pythia.readString("ExtraDimensionsLED:LambdaT = " + std::to_string(lambdaT)); //CMS observed lower limit 10TeV, so default 10000 GeV
    
    pythia.readString("PhaseSpace:pTHatMin = 20."); //increase this?
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13600.");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    pythia.init();

    pythia.settings.listAll();

    // Set up FastJet
    double R = 0.5;
    JetDefinition jet_def(antikt_algorithm, R);

    // Set up ROOT output
    std::string filename = "dijet_events_ADD_bsm_" + std::to_string(seed) + "_LT" + std::to_string(lambdaT) + ".root";
    TFile* outfile = new TFile(filename.c_str(), "RECREATE");
    TTree* tree = new TTree("dijets", "Dijet Tree");

    // Declare jet variables
    float jet1_px, jet1_py, jet1_pz, jet1_E, jet1_pt, jet1_eta, jet1_phi, jet1_mass;
    float jet2_px, jet2_py, jet2_pz, jet2_E, jet2_pt, jet2_eta, jet2_phi, jet2_mass;
    int jet1_partonID, jet2_partonID, jet1_nCons, jet2_nCons, jet1_nChargedCons, jet2_nChargedCons;
    
    // Declare jet constituent momentum, ID, and energy vectors 
    std::vector<float> jet1_Cons_Px, jet1_Cons_Py, jet1_Cons_Pz;
    std::vector<float> jet2_Cons_Px, jet2_Cons_Py, jet2_Cons_Pz;
    std::vector<int> jet1_Cons_ID;
    std::vector<int> jet2_Cons_ID;
    std::vector<float> jet1_Cons_E;
    std::vector<float> jet2_Cons_E;

    // Parton-level info (PDF inputs)
    int id1, id2;
    float x1, x2, Q2;

    // Create branches (same as before)
    tree->Branch("jet1_px", &jet1_px);
    tree->Branch("jet1_py", &jet1_py);
    tree->Branch("jet1_pz", &jet1_pz);
    tree->Branch("jet1_E", &jet1_E);
    tree->Branch("jet1_pt", &jet1_pt);
    tree->Branch("jet1_eta", &jet1_eta);
    tree->Branch("jet1_phi", &jet1_phi);
    tree->Branch("jet1_mass", &jet1_mass);
    tree->Branch("jet1_partonID", &jet1_partonID);
    tree->Branch("jet1_nCons", &jet1_nCons);
    tree->Branch("jet1_nChargedCons", &jet1_nChargedCons);
    tree->Branch("jet1_Cons_Px",&jet1_Cons_Px);
    tree->Branch("jet1_Cons_Py",&jet1_Cons_Py);
    tree->Branch("jet1_Cons_Pz",&jet1_Cons_Pz);
    tree->Branch("jet1_Cons_ID",&jet1_Cons_ID);
    tree->Branch("jet1_Cons_E", &jet1_Cons_E);


    tree->Branch("jet2_px", &jet2_px);
    tree->Branch("jet2_py", &jet2_py);
    tree->Branch("jet2_pz", &jet2_pz);
    tree->Branch("jet2_E", &jet2_E);
    tree->Branch("jet2_pt", &jet2_pt);
    tree->Branch("jet2_eta", &jet2_eta);
    tree->Branch("jet2_phi", &jet2_phi);
    tree->Branch("jet2_mass", &jet2_mass);
    tree->Branch("jet2_partonID", &jet2_partonID);
    tree->Branch("jet2_nCons", &jet2_nCons);
    tree->Branch("jet2_nChargedCons", &jet2_nChargedCons);
    tree->Branch("jet2_Cons_Px",&jet2_Cons_Px);
    tree->Branch("jet2_Cons_Py",&jet2_Cons_Py);
    tree->Branch("jet2_Cons_Pz",&jet2_Cons_Pz);
    tree->Branch("jet2_Cons_ID",&jet2_Cons_ID);
    tree->Branch("jet2_Cons_E", &jet2_Cons_E);

    tree->Branch("id1", &id1);
    tree->Branch("x1", &x1);
    tree->Branch("id2", &id2);
    tree->Branch("x2", &x2);
    tree->Branch("Q2", &Q2);

    int nSelected = 0;
    int nGenerated = 0;

    // Main event loop
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        nGenerated++;

        // Collect final visible particles and set user_index
        std::vector<PseudoJet> particles;
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];
            if (!p.isFinal() || !p.isVisible()) continue;
            PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
            pj.set_user_index(i);  // Store Pythia particle index
            particles.push_back(pj);
        }

        if (particles.size() < 2) continue;

        // Cluster jets
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(20.0));
        if (jets.size() < 2) continue;

        const PseudoJet& jet1 = jets[0];
        const PseudoJet& jet2 = jets[1];

        // Fiducial cuts
        if (jet1.pt() < 20.0 || jet2.pt() < 20.0) continue;
        if (jet1.eta() < 2.0 || jet1.eta() > 5.0) continue;
        if (jet2.eta() < 2.0 || jet2.eta() > 5.0) continue;

        nSelected++;

        // Fill jet variables
        jet1_px = jet1.px();
        jet1_py = jet1.py();
        jet1_pz = jet1.pz();
        jet1_E  = jet1.E();
        jet1_pt = jet1.pt();
        jet1_eta = jet1.eta();
        jet1_phi = jet1.phi_std();
        jet1_mass = jet1.m();
        jet1_partonID = getJetPartonID(jet1, pythia.event);
        jet1_nCons = jet1.constituents().size();

        // Count charged constituents in jet1
        jet1_nChargedCons = 0;
        for (const PseudoJet& c : jet1.constituents()) {
            int idx = c.user_index();
            if (idx < 0 || idx >= pythia.event.size()) continue;
            const Particle& p = pythia.event[idx];
            if (p.isCharged()) jet1_nChargedCons++;
        }

        jet2_px = jet2.px();
        jet2_py = jet2.py();
        jet2_pz = jet2.pz();
        jet2_E  = jet2.E();
        jet2_pt = jet2.pt();
        jet2_eta = jet2.eta();
        jet2_phi = jet2.phi_std();
        jet2_mass = jet2.m();
        jet2_partonID = getJetPartonID(jet2, pythia.event);
        jet2_nCons = jet2.constituents().size();

        // Count charged constituents in jet1
        jet2_nChargedCons = 0;
        for (const PseudoJet& c : jet2.constituents()) {
            int idx = c.user_index();
            if (idx < 0 || idx >= pythia.event.size()) continue;
            const Particle& p = pythia.event[idx];
            if (p.isCharged()) jet2_nChargedCons++;
        }

        // Fill PDF info
        id1 = pythia.info.id1pdf();
        id2 = pythia.info.id2pdf();
        x1  = pythia.info.x1pdf();
        x2  = pythia.info.x2pdf();
        Q2  = pythia.info.Q2Ren();

        // Clear vector values 
        jet1_Cons_Px.clear();
        jet1_Cons_Py.clear();
        jet1_Cons_Pz.clear();
        jet1_Cons_ID.clear();
        jet1_Cons_E.clear();
        jet2_Cons_Px.clear();
        jet2_Cons_Py.clear();        
        jet2_Cons_Pz.clear();        
        jet2_Cons_ID.clear();        
        jet2_Cons_E.clear();      
        
        // Fill jet constituent momenta and IDs 
        for (const PseudoJet& c : jet1.constituents()) {
            jet1_Cons_Px.emplace_back(c.px());
            jet1_Cons_Py.emplace_back(c.py());
            jet1_Cons_Pz.emplace_back(c.pz());
            jet1_Cons_E.push_back(c.E());
            int idx = c.user_index();
            if (idx >= 0 && idx < pythia.event.size()) {
                jet1_Cons_ID.push_back(pythia.event[idx].id());
            } else {
                jet1_Cons_ID.push_back(0);
            }
        }

        for (const PseudoJet& c : jet2.constituents()) {
            jet2_Cons_Px.emplace_back(c.px());
            jet2_Cons_Py.emplace_back(c.py());
            jet2_Cons_Pz.emplace_back(c.pz());
            jet2_Cons_E.push_back(c.E());
            int idx = c.user_index();
            if (idx >= 0 && idx < pythia.event.size()) {
                jet2_Cons_ID.push_back(pythia.event[idx].id());
            } else {
                jet2_Cons_ID.push_back(0);
            }
        }

        tree->Fill();
    }

    // Finalize
    pythia.stat();
    double sigmaTotal = pythia.info.sigmaGen(); 
    
    double sigmaDijet = sigmaTotal * (double)nSelected / (double)nGenerated;

    std::cout << "Dijet cross section estimate: " << sigmaDijet << " mb "  << std::endl; 
    
    outfile->Write();
    outfile->Close();

    return 0;
}
