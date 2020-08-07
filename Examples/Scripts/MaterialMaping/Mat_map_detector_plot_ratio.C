// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Mat_map_volume_plot.C"

/// Draw and save the ratio plots.

void plot_ratio(std::vector<TH2F*> Map_prop, std::vector<TH2F*> Map_geant, std::vector<int> detectors, const std::string& name){


  TH2F *Unit_Map_prop = (TH2F*) Map_prop[2]->Clone();
  Unit_Map_prop->Divide(Map_prop[2]);
  TH2F *Unit_Map_geant = (TH2F*) Map_geant[2]->Clone();
  Unit_Map_geant->Divide(Map_geant[2]);

  TH1D *Proj_eta_prop = (TH1D*) Map_prop[0]->ProjectionX()->Clone();
  Proj_eta_prop->Divide(Unit_Map_prop->ProjectionX());
  TH1D *Proj_eta_geant = (TH1D*) Map_geant[0]->ProjectionX()->Clone();
  Proj_eta_geant->Divide(Unit_Map_geant->ProjectionX());
    
  TH1D *Proj_phi_prop = (TH1D*) Map_prop[0]->ProjectionY()->Clone();
  Proj_phi_prop->Divide(Unit_Map_prop->ProjectionY());
  TH1D *Proj_phi_geant = (TH1D*) Map_geant[0]->ProjectionY()->Clone();
  Proj_phi_geant->Divide(Unit_Map_geant->ProjectionY());

  std::string sVol = "Detector volumes :";
  for(auto const& det: detectors) {
    sVol += " ";
    sVol += std::to_string(detectors);
  }
  TText *vol = new TText(.1, .95, sVol.c_str());
  vol->SetNDC();

  TCanvas *c1 = new TCanvas("c1","mat_X0",1200,1200);
  c1->SetRightMargin(0.14);
  c1->SetTopMargin(0.14);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);
  Map_prop[0]->Divide(Map_geant[0]);
  Map_prop[0]->GetZaxis()->SetTitle("X0 Val/Geant");
  Map_prop[0]->SetMaximum(2.);
  Map_prop[0]->Draw("COLZ");
  vol->Draw();
  c1->Print( (name+"ratio_X0.pdf").c_str());

  TCanvas *c2 = new TCanvas("c2", "mat_X0/eta", 1200, 1200);
  c2->SetRightMargin(0.14);
  c2->SetTopMargin(0.14);
  c2->SetLeftMargin(0.14);
  c2->SetBottomMargin(0.14);

  Proj_eta_prop->Divide(Proj_eta_geant);
  Proj_eta_prop->GetYaxis()->SetTitle("X0 Val/Geant");
  Proj_eta_prop->SetMarkerStyle(7);
  Proj_eta_prop->Draw("HIST PC");
  c2->Print((name + "ratio_X0_eta.pdf").c_str());

  TCanvas *c3 = new TCanvas("c3", "mat_X0/phi", 1200, 1200);
  c3->SetRightMargin(0.14);
  c3->SetTopMargin(0.14);
  c3->SetLeftMargin(0.14);
  c3->SetBottomMargin(0.14);
  Proj_phi_prop->Divide(Proj_phi_geant);
  Proj_phi_prop->GetYaxis()->SetTitle("X0 Val/Geant");
  Proj_phi_prop->SetMarkerStyle(7);
  Proj_phi_prop->Draw("HIST PC");
  c3->Print((name + "ratio_X0_phi.pdf").c_str());

  delete c1;
  delete c2;
  delete c3;
  delete vol;
  delete Unit_Map_prop;
  delete Unit_Map_geant;
} 


/// Plot the material ratio between the geantino scan and the map validation for each volume.
/// If a volume map json file is specify it is parse to associate name to the different volume id
/// nbprocess : number of parameter to be processed
/// name : name of the output directory.
/// name_prop : name of the output directory for the map valdation.
/// name_geant : name of the output directory for the geantino scan.
/// The map valdation and geantino scan plots are only saved if name_prop and name_geant are defined.
/// The parsing of the Json volume map file (use to associate the name to the volumes)
/// might not work with version of root newer that version 6.18.04

void Mat_map_volume_plot_ratio(std::string input_file_prop = "", std::string input_file_geant = "", std::vector<int> detectors, int nbprocess = -1, std::string name = "", std::string name_prop = "", std::string name_geant = ""){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::vector<TH2F*> volume_hist_prop;
  std::vector<TH2F*> volume_hist_geant;

  Fill(volume_hist_prop, input_file_prop, detectors, nbprocess);
  Fill(volume_hist_geant, input_file_geant, detectors, nbprocess);

  plot(volume_hist_prop, detectors, name_prop);
  plot(volume_hist_geant, detectors, name_geant);
  plot_ratio(volume_hist_prop, volume_hist_geant, detectors, name);
}
