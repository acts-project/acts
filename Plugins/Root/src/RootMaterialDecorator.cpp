// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootMaterialDecorator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TIterator.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/iter_find.hpp>

namespace Acts {
class ISurfaceMaterial;
class IVolumeMaterial;
}  // namespace Acts

using namespace Acts;

ActsPlugins::RootMaterialDecorator::RootMaterialDecorator(
    const ActsPlugins::RootMaterialDecorator::Config& config,
    Logging::Level level)
    : m_cfg(config),
      m_logger{getDefaultLogger("RootMaterialDecorator", level)} {
  // Validate the configuration
  if (m_cfg.accessorOptions.folderSurfaceNameBase.empty()) {
    throw std::invalid_argument("Missing surface folder name base");
  } else if (m_cfg.accessorOptions.folderVolumeNameBase.empty()) {
    throw std::invalid_argument("Missing volume folder name base");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  }

  // Set up ROOT I/O
  m_inputFile = TFile::Open(m_cfg.fileName.c_str());
  if (m_inputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName + "'");
  }

  ActsPlugins::RootMaterialMapIo accessor(m_cfg.accessorConfig,
                                          m_logger->clone("RootMaterialMapIo"));
  auto maps = accessor.read(*m_inputFile, m_cfg.accessorOptions);

  m_surfaceMaterialMap = std::move(maps.surfaceMaterials);

  // Get the list of keys from the file
  TList* tlist = m_inputFile->GetListOfKeys();
  auto tIter = tlist->MakeIterator();
  tIter->Reset();

  // Iterate over the keys in the file
  while (auto* key = static_cast<TKey*>(tIter->Next())) {
    // Remember the directory
    std::string tdName(key->GetName());

    std::vector<std::string> splitNames;
    iter_split(
        splitNames, tdName,
        boost::algorithm::first_finder(m_cfg.accessorConfig.volumePrefix));

    ACTS_VERBOSE("Processing directory: " << tdName);
    if (splitNames[0] == m_cfg.accessorOptions.folderVolumeNameBase) {
      // The volume material to be read in for this
      std::shared_ptr<const IVolumeMaterial> vMaterial = nullptr;
      // Volume key
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value volID = std::stoi(splitNames[0]);

      // Reconstruct the geometry ID
      auto geoID = GeometryIdentifier().withVolume(volID);
      ACTS_VERBOSE("GeometryIdentifier re-constructed as " << geoID);

      // Construct the names
      std::string nName = tdName + "/" + m_cfg.accessorConfig.nBinsHistName;
      std::string vName = tdName + "/" + m_cfg.accessorConfig.axisDirHistName;
      std::string oName =
          tdName + "/" + m_cfg.accessorConfig.axisBoundaryTypeHistName;
      std::string minName =
          tdName + "/" + m_cfg.accessorConfig.minRangeHistName;
      std::string maxName =
          tdName + "/" + m_cfg.accessorConfig.maxRangeHistName;
      std::string x0Name = tdName + "/" + m_cfg.accessorConfig.x0HistName;
      std::string l0Name = tdName + "/" + m_cfg.accessorConfig.l0HistName;
      std::string aName = tdName + "/" + m_cfg.accessorConfig.aHistName;
      std::string zName = tdName + "/" + m_cfg.accessorConfig.zHistName;
      std::string rhoName = tdName + "/" + m_cfg.accessorConfig.rhoHistName;

      // Get the histograms
      TH1F* n = dynamic_cast<TH1F*>(m_inputFile->Get(nName.c_str()));
      TH1F* v = dynamic_cast<TH1F*>(m_inputFile->Get(vName.c_str()));
      TH1F* o = dynamic_cast<TH1F*>(m_inputFile->Get(oName.c_str()));
      TH1F* min = dynamic_cast<TH1F*>(m_inputFile->Get(minName.c_str()));
      TH1F* max = dynamic_cast<TH1F*>(m_inputFile->Get(maxName.c_str()));
      TH1F* x0 = dynamic_cast<TH1F*>(m_inputFile->Get(x0Name.c_str()));
      TH1F* l0 = dynamic_cast<TH1F*>(m_inputFile->Get(l0Name.c_str()));
      TH1F* A = dynamic_cast<TH1F*>(m_inputFile->Get(aName.c_str()));
      TH1F* Z = dynamic_cast<TH1F*>(m_inputFile->Get(zName.c_str()));
      TH1F* rho = dynamic_cast<TH1F*>(m_inputFile->Get(rhoName.c_str()));

      // Only go on when you have all the material histograms
      if ((x0 != nullptr) && (l0 != nullptr) && (A != nullptr) &&
          (Z != nullptr) && (rho != nullptr)) {
        // If the bin information histograms are present the material is
        // either a 2D or a 3D grid
        if ((n != nullptr) && (v != nullptr) && (o != nullptr) &&
            (min != nullptr) && (max != nullptr)) {
          throw std::invalid_argument(
              "Interpolated volume material ROOT maps are no longer supported");
        } else {
          // Homogeneous material
          double dx0 = x0->GetBinContent(1);
          double dl0 = l0->GetBinContent(1);
          double da = A->GetBinContent(1);
          double dz = Z->GetBinContent(1);
          double drho = rho->GetBinContent(1);
          // Create material properties
          const auto material =
              Material::fromMassDensity(dx0, dl0, da, dz, drho);
          vMaterial = std::make_shared<HomogeneousVolumeMaterial>(material);
        }
      }
      ACTS_VERBOSE("Successfully read Material for : " << geoID);

      // Insert into the new collection
      m_volumeMaterialMap.try_emplace(geoID, std::move(vMaterial));
    }
  }
}

ActsPlugins::RootMaterialDecorator::~RootMaterialDecorator() {
  if (m_inputFile != nullptr) {
    m_inputFile->Close();
  }
}
