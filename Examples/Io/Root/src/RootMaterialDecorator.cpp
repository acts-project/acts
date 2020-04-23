// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootMaterialDecorator.hpp"

#include <Acts/Geometry/GeometryID.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Utilities/BinUtility.hpp>
#include <Acts/Utilities/BinningType.hpp>
#include <TFile.h>
#include <TH2F.h>
#include <TIterator.h>
#include <TKey.h>
#include <TList.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

FW::RootMaterialDecorator::RootMaterialDecorator(
    const FW::RootMaterialDecorator::Config& cfg)
    : m_cfg(cfg), m_inputFile(nullptr) {
  // Validate the configuration
  if (m_cfg.folderNameBase.empty()) {
    throw std::invalid_argument("Missing ROOT folder name");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  } else if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  } else if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing service name");
  }

  // Setup ROOT I/O
  m_inputFile = TFile::Open(m_cfg.fileName.c_str());
  if (!m_inputFile) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName);
  }

  // Get the list of keys from the file
  TList* tlist = m_inputFile->GetListOfKeys();
  auto tIter = tlist->MakeIterator();
  tIter->Reset();

  // Iterate over the keys in the file
  while (TKey* key = (TKey*)(tIter->Next())) {
    // The surface material to be read in for this
    std::shared_ptr<const Acts::ISurfaceMaterial> sMaterial = nullptr;

    // Remember the directory
    std::string tdName(key->GetName());

    ACTS_VERBOSE("Processing directory: " << tdName);

    // volume
    std::vector<std::string> splitNames;
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.voltag));
    boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
    Acts::GeometryID::Value volID = std::stoi(splitNames[0]);
    // boundary
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.boutag));
    boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
    Acts::GeometryID::Value bouID = std::stoi(splitNames[0]);
    // layer
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.laytag));
    boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
    Acts::GeometryID::Value layID = std::stoi(splitNames[0]);
    // approach
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.apptag));
    boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
    Acts::GeometryID::Value appID = std::stoi(splitNames[0]);
    // sensitive
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.sentag));
    Acts::GeometryID::Value senID = std::stoi(splitNames[1]);

    // Reconstruct the geometry ID
    Acts::GeometryID geoID;
    geoID.setVolume(volID);
    geoID.setBoundary(bouID);
    geoID.setLayer(layID);
    geoID.setApproach(appID);
    geoID.setSensitive(senID);
    ACTS_VERBOSE("GeometryID re-constructed as " << geoID);

    // Construct the names
    std::string nName = tdName + "/" + m_cfg.ntag;
    std::string vName = tdName + "/" + m_cfg.vtag;
    std::string oName = tdName + "/" + m_cfg.otag;
    std::string minName = tdName + "/" + m_cfg.mintag;
    std::string maxName = tdName + "/" + m_cfg.maxtag;
    std::string tName = tdName + "/" + m_cfg.ttag;
    std::string x0Name = tdName + "/" + m_cfg.x0tag;
    std::string l0Name = tdName + "/" + m_cfg.l0tag;
    std::string aName = tdName + "/" + m_cfg.atag;
    std::string zName = tdName + "/" + m_cfg.ztag;
    std::string rhoName = tdName + "/" + m_cfg.rhotag;

    // Get the histograms
    TH1F* n = dynamic_cast<TH1F*>(m_inputFile->Get(nName.c_str()));
    TH1F* v = dynamic_cast<TH1F*>(m_inputFile->Get(vName.c_str()));
    TH1F* o = dynamic_cast<TH1F*>(m_inputFile->Get(oName.c_str()));
    TH1F* min = dynamic_cast<TH1F*>(m_inputFile->Get(minName.c_str()));
    TH1F* max = dynamic_cast<TH1F*>(m_inputFile->Get(maxName.c_str()));
    TH2F* t = dynamic_cast<TH2F*>(m_inputFile->Get(tName.c_str()));
    TH2F* x0 = dynamic_cast<TH2F*>(m_inputFile->Get(x0Name.c_str()));
    TH2F* l0 = dynamic_cast<TH2F*>(m_inputFile->Get(l0Name.c_str()));
    TH2F* A = dynamic_cast<TH2F*>(m_inputFile->Get(aName.c_str()));
    TH2F* Z = dynamic_cast<TH2F*>(m_inputFile->Get(zName.c_str()));
    TH2F* rho = dynamic_cast<TH2F*>(m_inputFile->Get(rhoName.c_str()));

    // Only go on when you have all histograms
    if (n and v and o and min and max and t and x0 and l0 and A and Z and rho) {
      // Get the number of bins
      int nbins0 = t->GetNbinsX();
      int nbins1 = t->GetNbinsY();

      // The material matrix
      Acts::MaterialPropertiesMatrix materialMatrix(
          nbins1,
          Acts::MaterialPropertiesVector(nbins0, Acts::MaterialProperties()));

      // We need binned material properties
      if (nbins0 * nbins1 > 1) {
        // Fill the matrix first
        for (int ib0 = 1; ib0 <= nbins0; ++ib0) {
          for (int ib1 = 1; ib1 <= nbins1; ++ib1) {
            double dt = t->GetBinContent(ib0, ib1);
            if (dt > 0.) {
              double dx0 = x0->GetBinContent(ib0, ib1);
              double dl0 = l0->GetBinContent(ib0, ib1);
              double da = A->GetBinContent(ib0, ib1);
              double dz = Z->GetBinContent(ib0, ib1);
              double drho = rho->GetBinContent(ib0, ib1);
              // Create material properties
              materialMatrix[ib1 - 1][ib0 - 1] =
                  Acts::MaterialProperties(dx0, dl0, da, dz, drho, dt);
            }
          }
        }

        // Now reconstruct the bin untilities
        Acts::BinUtility bUtility;
        for (int ib = 1; ib < n->GetNbinsX() + 1; ++ib) {
          size_t nbins = size_t(n->GetBinContent(ib));
          Acts::BinningValue val = Acts::BinningValue(v->GetBinContent(ib));
          Acts::BinningOption opt = Acts::BinningOption(o->GetBinContent(ib));
          float rmin = min->GetBinContent(ib);
          float rmax = max->GetBinContent(ib);
          bUtility += Acts::BinUtility(nbins, rmin, rmax, opt, val);
        }
        ACTS_VERBOSE("Created " << bUtility);

        // Construct the binned material with the right bin utility
        sMaterial = std::make_shared<const Acts::BinnedSurfaceMaterial>(
            bUtility, std::move(materialMatrix));

      } else {
        // Only homogeneous material present
        double dt = t->GetBinContent(1, 1);
        double dx0 = x0->GetBinContent(1, 1);
        double dl0 = l0->GetBinContent(1, 1);
        double da = A->GetBinContent(1, 1);
        double dz = Z->GetBinContent(1, 1);
        double drho = rho->GetBinContent(1, 1);
        // Create and set the homogenous surface material
        sMaterial = std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            Acts::MaterialProperties(dx0, dl0, da, dz, drho, dt));
      }
    }
    ACTS_VERBOSE("Successfully read Material for : " << geoID);

    // Insert into the new collection
    m_surfaceMaterialMap.insert({geoID, std::move(sMaterial)});
  }
}

FW::RootMaterialDecorator::~RootMaterialDecorator() {
  m_inputFile->Close();
}
