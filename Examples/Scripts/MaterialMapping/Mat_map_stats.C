// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

void makeBinCountsDist(const char* filename = "material-maps.root")
{
  TFile* f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) return;

  TIter nextDir(f->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)nextDir())) {

    TObject* obj = key->ReadObj();

    if (!obj->InheritsFrom(TDirectory::Class())){
      continue;
    }

    TDirectory* dir = (TDirectory*)obj;

    TH2F* h2 = nullptr;
    dir->GetObject("binCounts", h2);

    if (!h2){
        std::cout<<"Histogram not found - continue"<<std::endl;
        continue;
    }

    std::cout << "\nProcessing: " << dir->GetName() << std::endl;

    TH1F hDist("binCountsDist",
               "BinCounts distribution;count;bins",
               100, 0, h2->GetMaximum() + 1);

    std::size_t nEmptyBins{0};

    for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= h2->GetNbinsY(); ++iy) {

        double val = h2->GetBinContent(ix, iy);
        if(val == 0){
            nEmptyBins+=1;
        }
        hDist.Fill(val);

      }
    }

    std::cout << "Mean = " << hDist.GetMean()
              << "  ,StdDev = " << hDist.GetStdDev()
              <<" ,Empty bins counted = "<< nEmptyBins
              << std::endl;

    TCanvas* c = new TCanvas();
    hDist.Draw();

    c->SaveAs(Form("%s_binCounts.png", dir->GetName()));
  }
}
