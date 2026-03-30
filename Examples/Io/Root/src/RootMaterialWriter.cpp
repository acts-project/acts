// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <ios>
#include <stdexcept>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

namespace ActsExamples {

RootMaterialWriter::RootMaterialWriter(const RootMaterialWriter::Config& config,
                                       Acts::Logging::Level level)
    : m_cfg(config),
      m_logger{Acts::getDefaultLogger("RootMaterialWriter", level)} {
  // Validate the configuration
  if (m_cfg.accessorOptions.folderSurfaceNameBase.empty()) {
    throw std::invalid_argument("Missing surface folder name base");
  } else if (m_cfg.accessorOptions.folderVolumeNameBase.empty()) {
    throw std::invalid_argument("Missing volume folder name base");
  } else if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file name");
  }

  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
}

void RootMaterialWriter::writeMaterial(
    const Acts::TrackingGeometryMaterial& detMaterial) {
  // Change to the output file
  m_outputFile->cd();

  const auto& [surfaceMaps, volumeMaps] = detMaterial;

  // Write the surface material maps
  ActsPlugins::RootMaterialMapIo accessor(m_cfg.accessorConfig,
                                          m_logger->clone("RootMaterialMapIo"));

  for (const auto& [geoId, sMap] : surfaceMaps) {
    // Get the Surface material
    accessor.write(*m_outputFile, geoId, *sMap, m_cfg.accessorOptions);
  }

  // Write the volume material maps
  for (auto& [key, value] : volumeMaps) {
    // Get the Volume material
    const Acts::IVolumeMaterial* vMaterial = value.get();
    if (vMaterial == nullptr) {
      ACTS_WARNING("No material for volume " << key << " skipping");
      continue;
    }

    // get the geometry ID
    Acts::GeometryIdentifier geoID = key;
    // decode the geometryID
    const auto gvolID = geoID.volume();

    // create the directory
    std::string tdName = m_cfg.accessorOptions.folderVolumeNameBase.c_str();
    tdName += m_cfg.accessorConfig.volumePrefix + std::to_string(gvolID);

    // create a new directory
    m_outputFile->mkdir(tdName.c_str());
    m_outputFile->cd(tdName.c_str());

    ACTS_VERBOSE("Writing out map at " << tdName);

    // understand what sort of material you have in mind
    auto bvMaterial3D = dynamic_cast<const Acts::InterpolatedMaterialMap<
        Acts::MaterialMapLookup<Acts::MaterialGrid3D>>*>(vMaterial);
    auto bvMaterial2D = dynamic_cast<const Acts::InterpolatedMaterialMap<
        Acts::MaterialMapLookup<Acts::MaterialGrid2D>>*>(vMaterial);

    int points = 1;
    if (bvMaterial3D != nullptr || bvMaterial2D != nullptr) {
      // Get the binning data
      std::vector<Acts::BinningData> binningData;
      if (bvMaterial3D != nullptr) {
        binningData = bvMaterial3D->binUtility().binningData();
        Acts::MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
        points = static_cast<int>(grid.size());
      } else {
        binningData = bvMaterial2D->binUtility().binningData();
        Acts::MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
        points = static_cast<int>(grid.size());
      }

      // 2-D or 3-D maps
      auto bins = static_cast<int>(binningData.size());
      auto fBins = static_cast<float>(bins);

      // The bin number information
      TH1F n(m_cfg.accessorConfig.nBinsHistName.c_str(), "bins; bin", bins,
             -0.5, fBins - 0.5);

      // The binning value information
      TH1F v(m_cfg.accessorConfig.axisDirHistName.c_str(),
             "binning values; bin", bins, -0.5, fBins - 0.5);

      // The binning option information
      TH1F o(m_cfg.accessorConfig.axisBoundaryTypeHistName.c_str(),
             "binning options; bin", bins, -0.5, fBins - 0.5);

      // The binning option information
      TH1F rmin(m_cfg.accessorConfig.minRangeHistName.c_str(), "min; bin", bins,
                -0.5, fBins - 0.5);

      // The binning option information
      TH1F rmax(m_cfg.accessorConfig.maxRangeHistName.c_str(), "max; bin", bins,
                -0.5, fBins - 0.5);

      // Now fill the histogram content
      for (const auto& [b, bData] : enumerate(binningData)) {
        // Fill: nbins, value, option, min, max
        n.SetBinContent(static_cast<int>(b),
                        static_cast<int>(binningData[b - 1].bins()));
        v.SetBinContent(static_cast<int>(b),
                        static_cast<int>(binningData[b - 1].binvalue));
        o.SetBinContent(static_cast<int>(b),
                        static_cast<int>(binningData[b - 1].option));
        rmin.SetBinContent(static_cast<int>(b), binningData[b - 1].min);
        rmax.SetBinContent(static_cast<int>(b), binningData[b - 1].max);
      }
      n.Write();
      v.Write();
      o.Write();
      rmin.Write();
      rmax.Write();
    }

    auto fPoints = static_cast<float>(points);
    TH1F x0(m_cfg.accessorConfig.x0HistName.c_str(), "X_{0} [mm] ;gridPoint",
            points, -0.5, fPoints - 0.5);
    TH1F l0(m_cfg.accessorConfig.l0HistName.c_str(),
            "#Lambda_{0} [mm] ;gridPoint", points, -0.5, fPoints - 0.5);
    TH1F A(m_cfg.accessorConfig.aHistName.c_str(), "X_{0} [mm] ;gridPoint",
           points, -0.5, fPoints - 0.5);
    TH1F Z(m_cfg.accessorConfig.zHistName.c_str(),
           "#Lambda_{0} [mm] ;gridPoint", points, -0.5, fPoints - 0.5);
    TH1F rho(m_cfg.accessorConfig.rhoHistName.c_str(),
             "#rho [g/mm^3] ;gridPoint", points, -0.5, fPoints - 0.5);
    // homogeneous volume
    if (points == 1) {
      auto mat = vMaterial->material({0, 0, 0});
      x0.SetBinContent(1, mat.X0());
      l0.SetBinContent(1, mat.L0());
      A.SetBinContent(1, mat.Ar());
      Z.SetBinContent(1, mat.Z());
      rho.SetBinContent(1, mat.massDensity());
    } else {
      // 3d grid volume
      if (bvMaterial3D != nullptr) {
        Acts::MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
        for (int point = 0; point < points; point++) {
          auto mat = Acts::Material(grid.at(point));
          if (!mat.isVacuum()) {
            x0.SetBinContent(point + 1, mat.X0());
            l0.SetBinContent(point + 1, mat.L0());
            A.SetBinContent(point + 1, mat.Ar());
            Z.SetBinContent(point + 1, mat.Z());
            rho.SetBinContent(point + 1, mat.massDensity());
          }
        }
      }
      // 2d grid volume
      else if (bvMaterial2D != nullptr) {
        Acts::MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
        for (int point = 0; point < points; point++) {
          auto mat = Acts::Material(grid.at(point));
          if (!mat.isVacuum()) {
            x0.SetBinContent(point + 1, mat.X0());
            l0.SetBinContent(point + 1, mat.L0());
            A.SetBinContent(point + 1, mat.Ar());
            Z.SetBinContent(point + 1, mat.Z());
            rho.SetBinContent(point + 1, mat.massDensity());
          }
        }
      }
    }
    x0.Write();
    l0.Write();
    A.Write();
    Z.Write();
    rho.Write();
  }
}

RootMaterialWriter::~RootMaterialWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

void RootMaterialWriter::write(const Acts::TrackingGeometry& tGeometry) {
  // Create a detector material map and loop recursively through it
  Acts::TrackingGeometryMaterial detMatMap;
  auto hVolume = tGeometry.highestTrackingVolume();
  if (hVolume != nullptr) {
    collectMaterial(*hVolume, detMatMap);
  }
  // Write the resulting map to the file
  writeMaterial(detMatMap);
}

void RootMaterialWriter::collectMaterial(
    const Acts::TrackingVolume& tVolume,
    Acts::TrackingGeometryMaterial& detMatMap) {
  // If the volume has volume material, write that
  if (tVolume.volumeMaterialPtr() != nullptr && m_cfg.processVolumes) {
    detMatMap.second[tVolume.geometryId()] = tVolume.volumeMaterialPtr();
  }

  // If confined layers exist, loop over them and collect the layer material
  if (tVolume.confinedLayers() != nullptr) {
    ACTS_VERBOSE("Collecting material for " << tVolume.volumeName()
                                            << " layers");
    for (auto& lay : tVolume.confinedLayers()->arrayObjects()) {
      collectMaterial(*lay, detMatMap);
    }
  }

  // If any of the boundary surfaces has material collect that
  if (m_cfg.processBoundaries) {
    for (auto& bou : tVolume.boundarySurfaces()) {
      const auto& bSurface = bou->surfaceRepresentation();
      if (bSurface.surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[bSurface.geometryId()] =
            bSurface.surfaceMaterialSharedPtr();
      }
    }
  }

  // If the volume has sub volumes, step down
  if (tVolume.confinedVolumes() != nullptr) {
    for (auto& tvol : tVolume.confinedVolumes()->arrayObjects()) {
      collectMaterial(*tvol, detMatMap);
    }
  }
}

void RootMaterialWriter::collectMaterial(
    const Acts::Layer& tLayer, Acts::TrackingGeometryMaterial& detMatMap) {
  // If the representing surface has material, collect it
  const auto& rSurface = tLayer.surfaceRepresentation();
  if (rSurface.surfaceMaterialSharedPtr() != nullptr &&
      m_cfg.processRepresenting) {
    detMatMap.first[rSurface.geometryId()] =
        rSurface.surfaceMaterialSharedPtr();
  }

  // Check the approach surfaces
  if (tLayer.approachDescriptor() != nullptr && m_cfg.processApproaches) {
    for (auto& aSurface : tLayer.approachDescriptor()->containedSurfaces()) {
      if (aSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[aSurface->geometryId()] =
            aSurface->surfaceMaterialSharedPtr();
      }
    }
  }

  // Check the sensitive surfaces
  if (tLayer.surfaceArray() != nullptr && m_cfg.processSensitives) {
    // sensitive surface loop
    for (auto& sSurface : tLayer.surfaceArray()->surfaces()) {
      if (sSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[sSurface->geometryId()] =
            sSurface->surfaceMaterialSharedPtr();
      }
    }
  }
}

}  // namespace ActsExamples
