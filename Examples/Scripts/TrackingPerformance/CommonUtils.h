// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <TColor.h>

// Helper function:
// function to set up the histogram style
template <typename hist_t>
void setHistStyle(hist_t* hist, short color = 1) {
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
  // hist->SetTitle("");
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
}

// Helper function:
// function to set up the efficiency histogram style
template <typename eff_t>
void setEffStyle(eff_t* eff, short color = 1) {
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(0.8);
  eff->SetLineWidth(2);
  eff->SetLineColor(color);
  eff->SetMarkerColor(color);
}

// Helper function:
// set color pallette
template <typename hist_t>
void adaptColorPalette(hist_t* h, float rmin, float rmax, float rgood,
                       float rwindow, int n) {
  // min - max is the range of the axis
  float rel_good = (rgood - rmin) / (rmax - rmin);
  float rel_window = rwindow / (rmax - rmin);

  // Stops are
  const int number = 5;
  double red[number] = {0., 0., 0., 1., 1.};
  double green[number] = {0., 1., 1., 1., 0.};
  double blue[number] = {1., 1., 0., 0., 0.};
  double stops[number] = {0., rel_good - rel_window, rel_good,
                          rel_good + rel_window, 1.};
  h->SetContour(n);

  TColor::CreateGradientColorTable(number, stops, red, green, blue, n);
}

// Helper function:
// increase eff range by a scale factor. Note that it assumes the eff has already been drawn 
template <typename eff_t>
void adaptEffRange(eff_t* eff, float minScale = 1, float maxScale = 1.1) {
  gPad->Update();
  auto ymin = gPad->GetUymin();
  auto ymax = gPad->GetUymax();
  auto graph = eff->GetPaintedGraph(); 
  graph->SetMinimum(ymin*minScale);
  graph->SetMaximum(ymax*maxScale); 
  gPad->Modified();
  gPad->Update();
}
