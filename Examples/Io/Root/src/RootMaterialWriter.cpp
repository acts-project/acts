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
#include "Acts/Material/GridVolumeMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
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

  const auto& volumeMaps = detMaterial.volumeMaterials;
  if (!detMaterial.keyedSurfaces.empty()) {
    throw std::invalid_argument(
        "ROOT material output does not support stable keys; use JSON");
  }

  // Write the surface material maps
  ActsPlugins::RootMaterialMapIo accessor(m_cfg.accessorConfig,
                                          m_logger->clone("RootMaterialMapIo"));

  accessor.write(*m_outputFile, detMaterial, m_cfg.accessorOptions);

  // Write the volume material maps
  for (auto& [key, value] : volumeMaps) {
    // Get the Volume material
    const Acts::IVolumeMaterial* vMaterial = value.get();
    if (vMaterial == nullptr) {
      ACTS_WARNING("No material for volume " << key << " skipping");
      continue;
    }
    if (dynamic_cast<const Acts::GridVolumeMaterial*>(vMaterial) != nullptr) {
      throw std::invalid_argument(
          "ROOT volume material output does not support grids; use JSON");
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

    const int points = 1;
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
    auto mat = vMaterial->material({0, 0, 0});
    x0.SetBinContent(1, mat.X0());
    l0.SetBinContent(1, mat.L0());
    A.SetBinContent(1, mat.Ar());
    Z.SetBinContent(1, mat.Z());
    rho.SetBinContent(1, mat.massDensity());
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
    detMatMap.volumeMaterials[tVolume.geometryId()] =
        tVolume.volumeMaterialPtr();
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
        detMatMap.surfaceMaterials[bSurface.geometryId()] =
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
    detMatMap.surfaceMaterials[rSurface.geometryId()] =
        rSurface.surfaceMaterialSharedPtr();
  }

  // Check the approach surfaces
  if (tLayer.approachDescriptor() != nullptr && m_cfg.processApproaches) {
    for (auto& aSurface : tLayer.approachDescriptor()->containedSurfaces()) {
      if (aSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.surfaceMaterials[aSurface->geometryId()] =
            aSurface->surfaceMaterialSharedPtr();
      }
    }
  }

  // Check the sensitive surfaces
  if (tLayer.surfaceArray() != nullptr && m_cfg.processSensitives) {
    // sensitive surface loop
    for (auto& sSurface : tLayer.surfaceArray()->surfaces()) {
      if (sSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.surfaceMaterials[sSurface->geometryId()] =
            sSurface->surfaceMaterialSharedPtr();
      }
    }
  }
}

}  // namespace ActsExamples
