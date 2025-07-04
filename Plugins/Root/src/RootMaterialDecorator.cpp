// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMaterialDecorator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <functional>
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

Acts::RootMaterialDecorator::RootMaterialDecorator(
    const Acts::RootMaterialDecorator::Config& config, Logging::Level level)
    : m_cfg(config),
      m_logger{getDefaultLogger("RootMaterialDecorator", level)} {
  // Validate the configuration
  if (m_cfg.folderSurfaceNameBase.empty()) {
    throw std::invalid_argument("Missing surface folder name base");
  } else if (m_cfg.folderVolumeNameBase.empty()) {
    throw std::invalid_argument("Missing volume folder name base");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  }

  // Setup ROOT I/O
  m_inputFile = TFile::Open(m_cfg.fileName.c_str());
  if (m_inputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName + "'");
  }

  // Get the list of keys from the file
  TList* tlist = m_inputFile->GetListOfKeys();
  auto tIter = tlist->MakeIterator();
  tIter->Reset();

  // Iterate over the keys in the file
  while (TKey* key = static_cast<TKey*>(tIter->Next())) {
    // Remember the directory
    std::string tdName(key->GetName());

    ACTS_VERBOSE("Processing directory: " << tdName);

    // volume
    std::vector<std::string> splitNames;
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.voltag));
    // Surface Material
    if (splitNames[0] == m_cfg.folderSurfaceNameBase) {
      // The surface material to be read in for this
      std::shared_ptr<const ISurfaceMaterial> sMaterial = nullptr;

      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value volID = std::stoi(splitNames[0]);
      // boundary
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.boutag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value bouID = std::stoi(splitNames[0]);
      // layer
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.laytag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value layID = std::stoi(splitNames[0]);
      // approach
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.apptag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value appID = std::stoi(splitNames[0]);
      // sensitive
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.sentag));
      GeometryIdentifier::Value senID = std::stoi(splitNames[1]);

      // Reconstruct the geometry ID
      auto geoID = GeometryIdentifier()
                       .withVolume(volID)
                       .withBoundary(bouID)
                       .withLayer(layID)
                       .withApproach(appID)
                       .withSensitive(senID);
      ACTS_VERBOSE("GeometryIdentifier re-constructed as " << geoID);

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
      auto n = dynamic_cast<TH1F*>(m_inputFile->Get(nName.c_str()));
      auto v = dynamic_cast<TH1F*>(m_inputFile->Get(vName.c_str()));
      auto o = dynamic_cast<TH1F*>(m_inputFile->Get(oName.c_str()));
      auto min = dynamic_cast<TH1F*>(m_inputFile->Get(minName.c_str()));
      auto max = dynamic_cast<TH1F*>(m_inputFile->Get(maxName.c_str()));
      auto t = dynamic_cast<TH2F*>(m_inputFile->Get(tName.c_str()));
      auto x0 = dynamic_cast<TH2F*>(m_inputFile->Get(x0Name.c_str()));
      auto l0 = dynamic_cast<TH2F*>(m_inputFile->Get(l0Name.c_str()));
      auto A = dynamic_cast<TH2F*>(m_inputFile->Get(aName.c_str()));
      auto Z = dynamic_cast<TH2F*>(m_inputFile->Get(zName.c_str()));
      auto rho = dynamic_cast<TH2F*>(m_inputFile->Get(rhoName.c_str()));

      std::vector<const TH1*> hists{n, v, o, min, max, t, x0, l0, A, Z, rho};

      // Only go on when you have all histograms
      if (std::ranges::all_of(
              hists, [](const auto* hist) { return hist != nullptr; })) {
        // Get the number of bins
        int nbins0 = t->GetNbinsX();
        int nbins1 = t->GetNbinsY();

        // The material matrix
        MaterialSlabMatrix materialMatrix(
            nbins1, MaterialSlabVector(nbins0, MaterialSlab::Nothing()));

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
                const auto material =
                    Material::fromMassDensity(dx0, dl0, da, dz, drho);
                materialMatrix[ib1 - 1][ib0 - 1] = MaterialSlab(material, dt);
              }
            }
          }

          // Now reconstruct the bin untilities
          BinUtility bUtility;
          for (int ib = 1; ib < n->GetNbinsX() + 1; ++ib) {
            std::size_t nbins = static_cast<std::size_t>(n->GetBinContent(ib));
            auto val = static_cast<AxisDirection>(v->GetBinContent(ib));
            auto opt = static_cast<BinningOption>(o->GetBinContent(ib));
            float rmin = static_cast<float>(min->GetBinContent(ib));
            float rmax = static_cast<float>(max->GetBinContent(ib));
            bUtility += BinUtility(nbins, rmin, rmax, opt, val);
          }
          ACTS_VERBOSE("Created " << bUtility);

          // Construct the binned material with the right bin utility
          sMaterial = std::make_shared<const BinnedSurfaceMaterial>(
              bUtility, std::move(materialMatrix));

        } else {
          // Only homogeneous material present
          double dt = t->GetBinContent(1, 1);
          double dx0 = x0->GetBinContent(1, 1);
          double dl0 = l0->GetBinContent(1, 1);
          double da = A->GetBinContent(1, 1);
          double dz = Z->GetBinContent(1, 1);
          double drho = rho->GetBinContent(1, 1);
          // Create and set the homogeneous surface material
          const auto material =
              Material::fromMassDensity(dx0, dl0, da, dz, drho);
          sMaterial = std::make_shared<const HomogeneousSurfaceMaterial>(
              MaterialSlab(material, dt));
        }
      }
      ACTS_VERBOSE("Successfully read Material for : " << geoID);

      // Insert into the new collection
      m_surfaceMaterialMap.insert({geoID, std::move(sMaterial)});

    } else if (splitNames[0] == m_cfg.folderVolumeNameBase) {
      // The volume material to be read in for this
      std::shared_ptr<const IVolumeMaterial> vMaterial = nullptr;
      // Volume key
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value volID = std::stoi(splitNames[0]);

      // Reconstruct the geometry ID
      auto geoID = GeometryIdentifier().withVolume(volID);
      ACTS_VERBOSE("GeometryIdentifier re-constructed as " << geoID);

      // Construct the names
      std::string nName = tdName + "/" + m_cfg.ntag;
      std::string vName = tdName + "/" + m_cfg.vtag;
      std::string oName = tdName + "/" + m_cfg.otag;
      std::string minName = tdName + "/" + m_cfg.mintag;
      std::string maxName = tdName + "/" + m_cfg.maxtag;
      std::string x0Name = tdName + "/" + m_cfg.x0tag;
      std::string l0Name = tdName + "/" + m_cfg.l0tag;
      std::string aName = tdName + "/" + m_cfg.atag;
      std::string zName = tdName + "/" + m_cfg.ztag;
      std::string rhoName = tdName + "/" + m_cfg.rhotag;

      // Get the histograms
      auto n = dynamic_cast<TH1F*>(m_inputFile->Get(nName.c_str()));
      auto v = dynamic_cast<TH1F*>(m_inputFile->Get(vName.c_str()));
      auto o = dynamic_cast<TH1F*>(m_inputFile->Get(oName.c_str()));
      auto min = dynamic_cast<TH1F*>(m_inputFile->Get(minName.c_str()));
      auto max = dynamic_cast<TH1F*>(m_inputFile->Get(maxName.c_str()));
      auto x0 = dynamic_cast<TH1F*>(m_inputFile->Get(x0Name.c_str()));
      auto l0 = dynamic_cast<TH1F*>(m_inputFile->Get(l0Name.c_str()));
      auto A = dynamic_cast<TH1F*>(m_inputFile->Get(aName.c_str()));
      auto Z = dynamic_cast<TH1F*>(m_inputFile->Get(zName.c_str()));
      auto rho = dynamic_cast<TH1F*>(m_inputFile->Get(rhoName.c_str()));

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
          BinUtility bUtility;
          for (int ib = 1; ib < dim + 1; ++ib) {
            std::size_t nbins = static_cast<std::size_t>(n->GetBinContent(ib));
            auto val = static_cast<AxisDirection>(v->GetBinContent(ib));
            auto opt = static_cast<BinningOption>(o->GetBinContent(ib));
            float rmin = min->GetBinContent(ib);
            float rmax = max->GetBinContent(ib);
            bUtility += BinUtility(nbins, rmin, rmax, opt, val);
          }
          ACTS_VERBOSE("Created " << bUtility);

          if (dim == 2) {
            // 2D Grid material
            std::function<Vector2(Vector3)> transfoGlobalToLocal;
            Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

            Grid2D::point_t gMin = grid.minPosition();
            Grid2D::point_t gMax = grid.maxPosition();
            Grid2D::index_t gNBins = grid.numLocalBins();

            EAxis axis1(gMin[0], gMax[0], gNBins[0]);
            EAxis axis2(gMin[1], gMax[1], gNBins[1]);

            // Build the grid and fill it with data
            MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

            for (int p = 1; p <= points; p++) {
              double dx0 = x0->GetBinContent(p);
              double dl0 = l0->GetBinContent(p);
              double da = A->GetBinContent(p);
              double dz = Z->GetBinContent(p);
              double drho = rho->GetBinContent(p);
              // Create material properties
              const auto material =
                  Material::fromMassDensity(dx0, dl0, da, dz, drho);
              mGrid.at(p - 1) = material.parameters();
            }
            MaterialMapper<MaterialGrid2D> matMap(transfoGlobalToLocal, mGrid);
            vMaterial = std::make_shared<
                InterpolatedMaterialMap<MaterialMapper<MaterialGrid2D>>>(
                std::move(matMap), bUtility);
          } else if (dim == 3) {
            // 3D Grid material
            std::function<Vector3(Vector3)> transfoGlobalToLocal;
            Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

            Grid3D::point_t gMin = grid.minPosition();
            Grid3D::point_t gMax = grid.maxPosition();
            Grid3D::index_t gNBins = grid.numLocalBins();

            EAxis axis1(gMin[0], gMax[0], gNBins[0]);
            EAxis axis2(gMin[1], gMax[1], gNBins[1]);
            EAxis axis3(gMin[2], gMax[2], gNBins[2]);

            // Build the grid and fill it with data
            MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

            for (int p = 1; p <= points; p++) {
              double dx0 = x0->GetBinContent(p);
              double dl0 = l0->GetBinContent(p);
              double da = A->GetBinContent(p);
              double dz = Z->GetBinContent(p);
              double drho = rho->GetBinContent(p);
              // Create material properties
              const auto material = Material::fromMassDensity(
                  static_cast<float>(dx0), static_cast<float>(dl0),
                  static_cast<float>(da), static_cast<float>(dz),
                  static_cast<float>(drho));
              mGrid.at(p - 1) = material.parameters();
            }
            MaterialMapper<MaterialGrid3D> matMap(transfoGlobalToLocal, mGrid);
            vMaterial = std::make_shared<
                InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>>(
                std::move(matMap), bUtility);
          }
        } else {
          // Homogeneous material
          double dx0 = x0->GetBinContent(1);
          double dl0 = l0->GetBinContent(1);
          double da = A->GetBinContent(1);
          double dz = Z->GetBinContent(1);
          double drho = rho->GetBinContent(1);
          // Create material properties
          const auto material = Material::fromMassDensity(
              static_cast<float>(dx0), static_cast<float>(dl0),
              static_cast<float>(da), static_cast<float>(dz),
              static_cast<float>(drho));
          vMaterial = std::make_shared<HomogeneousVolumeMaterial>(material);
        }
      }
      ACTS_VERBOSE("Successfully read Material for : " << geoID);

      // Insert into the new collection
      m_volumeMaterialMap.insert({geoID, std::move(vMaterial)});

      // Incorrect FolderName value
    } else {
      ACTS_ERROR(
          "Invalid FolderName, does not match any entry in the root file");
    }
  }
}

Acts::RootMaterialDecorator::~RootMaterialDecorator() {
  if (m_inputFile != nullptr) {
    m_inputFile->Close();
  }
}
