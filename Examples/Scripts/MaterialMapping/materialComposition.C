// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <map>
#include <string>
#include <tuple>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TTree.h>

struct MaterialHistograms {
  TProfile* x0_vs_eta = nullptr;
  TProfile* l0_vs_eta = nullptr;

  TProfile* x0_vs_phi = nullptr;
  TProfile* l0_vs_phi = nullptr;

  float s_x0 = 0.;
  float s_l0 = 0.;

  /// Material histogram constructor
  ///
  /// @param iA the atomic number
  MaterialHistograms(const std::string& name, unsigned int iA,
                     unsigned int bins, float eta) {
    std::string x0NameEta =
        (iA == 0) ? name + std::string("_x0_vs_eta_all")
                  : name + std::string("_x0_vs_eta_A") + std::to_string(iA);
    std::string l0NameEta =
        (iA == 0) ? name + std::string("_l0_vs_eta_all")
                  : name + std::string("_l0_vs_eta_A") + std::to_string(iA);

    x0_vs_eta = new TProfile(x0NameEta.c_str(), "X_{0} vs. #eta", bins, -eta,
                             eta, 0., 5.);
    l0_vs_eta = new TProfile(l0NameEta.c_str(), "L_{0} vs. #eta", bins, -eta,
                             eta, 0., 5.);

    std::string x0NamePhi =
        (iA == 0) ? name + std::string("_x0_vs_phi_all")
                  : name + std::string("_x0_vs_phi_A") + std::to_string(iA);
    std::string l0NamePhi =
        (iA == 0) ? name + std::string("_l0_vs_phi_all")
                  : name + std::string("_l0_vs_phi_A") + std::to_string(iA);

    x0_vs_phi = new TProfile(x0NamePhi.c_str(), "X_{0} vs. #phi", bins, -M_PI,
                             M_PI, 0., 5.);
    l0_vs_phi = new TProfile(l0NamePhi.c_str(), "L_{0} vs. #phi", bins, -M_PI,
                             M_PI, 0., 5.);
  }

  MaterialHistograms() = default;

  /// This fills the event into the histograms
  /// and clears the cache accordingly
  ///
  /// @param eta the pseudorapidity value
  /// @param phi the phi value
  ///
  void fillAndClear(float eta, float phi) {
    
    x0_vs_eta->Fill(eta, s_x0);
    l0_vs_eta->Fill(eta, s_l0);

    x0_vs_phi->Fill(phi, s_x0);
    l0_vs_phi->Fill(phi, s_l0);

    s_x0 = 0.;
    s_l0 = 0.;
  }

  /// Write out the histograms, the TDirectory needs
  /// to be set before
  ///
  /// Histrograms with no contribution will not be 
  /// written to file.
  void write() {
    if (x0_vs_eta->GetMaximum() > 0.) {
      x0_vs_eta->Write();
      l0_vs_eta->Write();

      x0_vs_phi->Write();
      l0_vs_phi->Write();
    }
  }
};

using Region = std::tuple<std::string, float, float, float, float>;

/// Plot the material composition
///
/// @param inFile the input root file
/// @param treeNAme the input tree name (default: 'trackstates)
/// @param outFile the output root file
/// @param bins the number of bins
/// @param eta the eta range
void materialComposition(const std::string& inFile, const std::string& treeName,
                         const std::string& outFile, unsigned int bins,
                         float eta, const std::vector<Region>& regions) {
  // Open the input file & get the tree
  auto inputFile = TFile::Open(inFile.c_str());
  auto inputTree = dynamic_cast<TTree*>(inputFile->Get(treeName.c_str()));
  if (inputTree != nullptr) {
    // Get the different atomic numbers
    TCanvas* materialCanvas =
        new TCanvas("materialCanvas", "Materials", 100, 100, 620, 400);
    materialCanvas->cd();
    // Draw all the atomic elements & get the histogram
    inputTree->Draw("mat_A>>hA(100,0.5,100.5)");
    TH1F* histA = dynamic_cast<TH1F*>(gDirectory->Get("hA"));
    histA->Draw();

    auto outputFile = TFile::Open(outFile.c_str(), "recreate");

    float v_eta;
    float v_phi;
    std::vector<float>* stepLength = new std::vector<float>;
    std::vector<float>* stepX0 = new std::vector<float>;
    std::vector<float>* stepL0 = new std::vector<float>;
    std::vector<float>* stepA = new std::vector<float>;

    std::vector<float>* startX = new std::vector<float>;
    std::vector<float>* startY = new std::vector<float>;
    std::vector<float>* startZ = new std::vector<float>;

    std::vector<float>* endX = new std::vector<float>;
    std::vector<float>* endY = new std::vector<float>;
    std::vector<float>* endZ = new std::vector<float>;

    inputTree->SetBranchAddress("v_eta", &v_eta);
    inputTree->SetBranchAddress("v_phi", &v_phi);
    inputTree->SetBranchAddress("mat_step_length", &stepLength);

    inputTree->SetBranchAddress("mat_X0", &stepX0);
    inputTree->SetBranchAddress("mat_L0", &stepL0);
    inputTree->SetBranchAddress("mat_A", &stepA);

    inputTree->SetBranchAddress("mat_sx", &startX);
    inputTree->SetBranchAddress("mat_sy", &startY);
    inputTree->SetBranchAddress("mat_sz", &startZ);

    inputTree->SetBranchAddress("mat_ex", &endX);
    inputTree->SetBranchAddress("mat_ey", &endY);
    inputTree->SetBranchAddress("mat_ez", &endZ);

    // Loop over all entries ---------------
    unsigned int entries = inputTree->GetEntries();

#ifdef BOOST_AVAILABLE
    std::cout << "*** Event Loop: " << std::endl;
    progress_display event_loop_progress(entries * regions.size());
#endif

    // Loop of the regions
    for (auto& r : regions) {
      std::string rName = std::get<0>(r);

      // The material histograms ordered by atomic mass
      std::map<unsigned int, MaterialHistograms> mCache;

      // The material histograms for all
      mCache[0] = MaterialHistograms(rName, 0, bins, eta);
      for (unsigned int ib = 1; ib <= 100; ++ib) {
        if (histA->GetBinContent(ib) > 0.) {
          mCache[ib] = MaterialHistograms(rName, ib, bins, eta);
        }
      }

      for (unsigned int ie = 0; ie < entries; ++ie) {
        // Get the entry
        inputTree->GetEntry(ie);

#ifdef BOOST_AVAILABLE
        ++event_loop_progress;
#endif

        // Accumulate the material per track
        size_t steps = stepLength->size();
        for (unsigned int is = 0; is < steps; ++is) {
          float sX = startX->at(is);
          float sY = startY->at(is);
          float sZ = startZ->at(is);
          float sR = sqrt(sX * sX + sY * sY);

          float eX = endX->at(is);
          float eY = endY->at(is);
          float eZ = endZ->at(is);
          float eR = sqrt(eX * eX + eY * eY);

          float minR = std::get<1>(r);
          float maxR = std::get<2>(r);
          float minZ = std::get<3>(r);
          float maxZ = std::get<4>(r);

          if (minR > sR or minZ > sZ or maxR < eR or maxZ < eZ) {
            continue;
          }

          float step = stepLength->at(is);
          float X0 = stepX0->at(is);
          float L0 = stepL0->at(is);

          // The integral one
          auto& all = mCache[0];
          all.s_x0 += step / X0;
          all.s_l0 += step / L0;

          unsigned int sA = histA->FindBin(stepA->at(is));
          // The current one
          auto& current = mCache[sA];
          current.s_x0 += step / X0;
          current.s_l0 += step / L0;
        }
        // Fill the profiles and clear the cache
        for (auto& [key, cache] : mCache) {
          cache.fillAndClear(v_eta, v_phi);
        }
      }
      // -----------------------------------

      for (auto [key, cache] : mCache) {
        cache.write();
      }
    }
    outputFile->Close();
  }
}
