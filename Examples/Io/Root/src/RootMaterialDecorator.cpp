// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousVolumeMaterial.hpp>
#include <Acts/Utilities/AxisDefinitions.hpp>
#include <Acts/Utilities/BinUtility.hpp>
#include <Acts/Utilities/BinningType.hpp>

#include <algorithm>
#include <cstdio>
#include <functional>
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

ActsExamples::RootMaterialDecorator::RootMaterialDecorator(
    const ActsExamples::RootMaterialDecorator::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      m_logger{Acts::getDefaultLogger("RootMaterialDecorator", level)} {
  // Validate the configuration
  if (m_cfg.accessorOptions.folderSurfaceNameBase.empty()) {
    throw std::invalid_argument("Missing surface folder name base");
  } else if (m_cfg.accessorOptions.folderVolumeNameBase.empty()) {
    throw std::invalid_argument("Missing volume folder name base");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  }

  // Setup ROOT I/O
  m_inputFile = TFile::Open(m_cfg.fileName.c_str());
  if (m_inputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName + "'");
  }

  Acts::RootMaterialMapIo accessor(m_cfg.accessorConfig,
                                   m_logger->clone("RootMaterialMapIo"));
  auto [surfaceMaps, volumeMaps] =
      accessor.read(*m_inputFile, m_cfg.accessorOptions);

  m_surfaceMaterialMap = std::move(surfaceMaps);

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
      std::shared_ptr<const Acts::IVolumeMaterial> vMaterial = nullptr;
      // Volume key
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      Acts::GeometryIdentifier::Value volID = std::stoi(splitNames[0]);

      // Reconstruct the geometry ID
      auto geoID = Acts::GeometryIdentifier().withVolume(volID);
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
        // Get the number of grid points
        int points = x0->GetNbinsX();
        // If the bin information histograms are present the material is
        // either a 2D or a 3D grid
        if ((n != nullptr) && (v != nullptr) && (o != nullptr) &&
            (min != nullptr) && (max != nullptr)) {
          // Dimension of the grid
          int dim = n->GetNbinsX();
          // Now reconstruct the bin untilities
          Acts::BinUtility bUtility;
          for (int ib = 1; ib < dim + 1; ++ib) {
            std::size_t nbins = static_cast<std::size_t>(n->GetBinContent(ib));
            Acts::AxisDirection val =
                static_cast<Acts::AxisDirection>(v->GetBinContent(ib));
            Acts::BinningOption opt =
                static_cast<Acts::BinningOption>(o->GetBinContent(ib));
            float rmin = min->GetBinContent(ib);
            float rmax = max->GetBinContent(ib);
            bUtility += Acts::BinUtility(nbins, rmin, rmax, opt, val);
          }
          ACTS_VERBOSE("Created " << bUtility);

          if (dim == 2) {
            // 2D Grid material
            std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
            Acts::Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

            Acts::Grid2D::point_t gMin = grid.minPosition();
            Acts::Grid2D::point_t gMax = grid.maxPosition();
            Acts::Grid2D::index_t gNBins = grid.numLocalBins();

            Acts::EAxis axis1(gMin[0], gMax[0], gNBins[0]);
            Acts::EAxis axis2(gMin[1], gMax[1], gNBins[1]);

            // Build the grid and fill it with data
            Acts::MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

            for (int p = 1; p <= points; p++) {
              double dx0 = x0->GetBinContent(p);
              double dl0 = l0->GetBinContent(p);
              double da = A->GetBinContent(p);
              double dz = Z->GetBinContent(p);
              double drho = rho->GetBinContent(p);
              // Create material properties
              const auto material =
                  Acts::Material::fromMassDensity(dx0, dl0, da, dz, drho);
              mGrid.at(p - 1) = material.parameters();
            }
            Acts::MaterialMapper<Acts::MaterialGrid2D> matMap(
                transfoGlobalToLocal, mGrid);
            vMaterial = std::make_shared<Acts::InterpolatedMaterialMap<
                Acts::MaterialMapper<Acts::MaterialGrid2D>>>(std::move(matMap),
                                                             bUtility);
          } else if (dim == 3) {
            // 3D Grid material
            std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
            Acts::Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

            Acts::Grid3D::point_t gMin = grid.minPosition();
            Acts::Grid3D::point_t gMax = grid.maxPosition();
            Acts::Grid3D::index_t gNBins = grid.numLocalBins();

            Acts::EAxis axis1(gMin[0], gMax[0], gNBins[0]);
            Acts::EAxis axis2(gMin[1], gMax[1], gNBins[1]);
            Acts::EAxis axis3(gMin[2], gMax[2], gNBins[2]);

            // Build the grid and fill it with data
            Acts::MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

            for (int p = 1; p <= points; p++) {
              double dx0 = x0->GetBinContent(p);
              double dl0 = l0->GetBinContent(p);
              double da = A->GetBinContent(p);
              double dz = Z->GetBinContent(p);
              double drho = rho->GetBinContent(p);
              // Create material properties
              const auto material =
                  Acts::Material::fromMassDensity(dx0, dl0, da, dz, drho);
              mGrid.at(p - 1) = material.parameters();
            }
            Acts::MaterialMapper<Acts::MaterialGrid3D> matMap(
                transfoGlobalToLocal, mGrid);
            vMaterial = std::make_shared<Acts::InterpolatedMaterialMap<
                Acts::MaterialMapper<Acts::MaterialGrid3D>>>(std::move(matMap),
                                                             bUtility);
          }
        } else {
          // Homogeneous material
          double dx0 = x0->GetBinContent(1);
          double dl0 = l0->GetBinContent(1);
          double da = A->GetBinContent(1);
          double dz = Z->GetBinContent(1);
          double drho = rho->GetBinContent(1);
          // Create material properties
          const auto material =
              Acts::Material::fromMassDensity(dx0, dl0, da, dz, drho);
          vMaterial =
              std::make_shared<Acts::HomogeneousVolumeMaterial>(material);
        }
      }
      ACTS_VERBOSE("Successfully read Material for : " << geoID);

      // Insert into the new collection
      m_volumeMaterialMap.insert({geoID, std::move(vMaterial)});
    }
  }
}

ActsExamples::RootMaterialDecorator::~RootMaterialDecorator() {
  if (m_inputFile != nullptr) {
    m_inputFile->Close();
  }
}
