// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootMaterialWriter.hpp"

#include <Acts/Geometry/GeometryID.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <TFile.h>
#include <TH2F.h>
#include <ios>
#include <iostream>
#include <stdexcept>

FW::RootMaterialWriter::RootMaterialWriter(
    const FW::RootMaterialWriter::Config& cfg)
    : m_cfg(cfg) {
  // Validate the configuration
  if (m_cfg.folderNameBase.empty()) {
    throw std::invalid_argument("Missing folder name base");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  } else if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  } else if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing service name");
  }
}

void FW::RootMaterialWriter::write(
    const Acts::DetectorMaterialMaps& detMaterial) {
  // Setup ROOT I/O
  TFile* outputFile = TFile::Open(m_cfg.fileName.c_str(), "recreate");
  if (!outputFile) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName);
  }

  // Change to the output file
  outputFile->cd();

  auto& surfaceMaps = detMaterial.first;
  for (auto& [key, value] : surfaceMaps) {
    // Get the Surface material
    const Acts::ISurfaceMaterial* sMaterial = value.get();

    // get the geometry ID
    Acts::GeometryID geoID = key;
    // decode the geometryID
    const auto gvolID = geoID.volume();
    const auto gbouID = geoID.boundary();
    const auto glayID = geoID.layer();
    const auto gappID = geoID.approach();
    const auto gsenID = geoID.sensitive();
    // create the directory
    std::string tdName = m_cfg.folderNameBase.c_str();
    tdName += m_cfg.voltag + std::to_string(gvolID);
    tdName += m_cfg.boutag + std::to_string(gbouID);
    tdName += m_cfg.laytag + std::to_string(glayID);
    tdName += m_cfg.apptag + std::to_string(gappID);
    tdName += m_cfg.sentag + std::to_string(gsenID);
    // create a new directory
    outputFile->mkdir(tdName.c_str());
    outputFile->cd(tdName.c_str());

    ACTS_VERBOSE("Writing out map at " << tdName);

    size_t bins0 = 1, bins1 = 1;
    // understand what sort of material you have in mind
    const Acts::BinnedSurfaceMaterial* bsm =
        dynamic_cast<const Acts::BinnedSurfaceMaterial*>(sMaterial);
    if (bsm) {
      // overwrite the bin numbers
      bins0 = bsm->binUtility().bins(0);
      bins1 = bsm->binUtility().bins(1);

      // Get the binning data
      auto& binningData = bsm->binUtility().binningData();
      // 1-D or 2-D maps
      size_t binningBins = binningData.size();

      // The bin number information
      TH1F* n = new TH1F(m_cfg.ntag.c_str(), "bins; bin", binningBins, -0.5,
                         binningBins - 0.5);

      // The binning value information
      TH1F* v = new TH1F(m_cfg.vtag.c_str(), "binning values; bin", binningBins,
                         -0.5, binningBins - 0.5);

      // The binning option information
      TH1F* o = new TH1F(m_cfg.otag.c_str(), "binning options; bin",
                         binningBins, -0.5, binningBins - 0.5);

      // The binning option information
      TH1F* min = new TH1F(m_cfg.mintag.c_str(), "min; bin", binningBins, -0.5,
                           binningBins - 0.5);

      // The binning option information
      TH1F* max = new TH1F(m_cfg.maxtag.c_str(), "max; bin", binningBins, -0.5,
                           binningBins - 0.5);

      // Now fill the histogram content
      size_t b = 1;
      for (auto bData : binningData) {
        // Fill: nbins, value, option, min, max
        n->SetBinContent(b, int(binningData[b - 1].bins()));
        v->SetBinContent(b, int(binningData[b - 1].binvalue));
        o->SetBinContent(b, int(binningData[b - 1].option));
        min->SetBinContent(b, binningData[b - 1].min);
        max->SetBinContent(b, binningData[b - 1].max);
        ++b;
      }
      n->Write();
      v->Write();
      o->Write();
      min->Write();
      max->Write();
    }

    TH2F* t = new TH2F(m_cfg.ttag.c_str(), "thickness [mm] ;b0 ;b1", bins0,
                       -0.5, bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
    TH2F* x0 = new TH2F(m_cfg.x0tag.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5,
                        bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
    TH2F* l0 = new TH2F(m_cfg.l0tag.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0,
                        -0.5, bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
    TH2F* A = new TH2F(m_cfg.atag.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5,
                       bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
    TH2F* Z = new TH2F(m_cfg.ztag.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0,
                       -0.5, bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
    TH2F* rho = new TH2F(m_cfg.rhotag.c_str(), "#rho [g/mm^3] ;b0 ;b1", bins0,
                         -0.5, bins0 - 0.5, bins1, -0.5, bins1 - 0.5);

    // loop over the material and fill
    for (size_t b0 = 0; b0 < bins0; ++b0) {
      for (size_t b1 = 0; b1 < bins1; ++b1) {
        // get the material for the bin
        auto& mat = sMaterial->materialProperties(b0, b1);
        if (mat) {
          t->SetBinContent(b0 + 1, b1 + 1, mat.thickness());
          x0->SetBinContent(b0 + 1, b1 + 1, mat.material().X0());
          l0->SetBinContent(b0 + 1, b1 + 1, mat.material().L0());
          A->SetBinContent(b0 + 1, b1 + 1, mat.material().Ar());
          Z->SetBinContent(b0 + 1, b1 + 1, mat.material().Z());
          rho->SetBinContent(b0 + 1, b1 + 1, mat.material().massDensity());
        }
      }
    }
    t->Write();
    x0->Write();
    l0->Write();
    A->Write();
    Z->Write();
    rho->Write();
  }

  outputFile->Close();
}

void FW::RootMaterialWriter::write(const Acts::TrackingGeometry& tGeometry) {
  // Create a detector material map and loop recursively through it
  Acts::DetectorMaterialMaps detMatMap;
  auto hVolume = tGeometry.highestTrackingVolume();
  if (hVolume != nullptr) {
    collectMaterial(*hVolume, detMatMap);
  }
  // Write the resulting map to the file
  write(detMatMap);
}

void FW::RootMaterialWriter::collectMaterial(
    const Acts::TrackingVolume& tVolume,
    Acts::DetectorMaterialMaps& detMatMap) {
  // If the volume has volume material, write that
  if (tVolume.volumeMaterialSharedPtr() != nullptr and m_cfg.processVolumes) {
    detMatMap.second[tVolume.geoID()] = tVolume.volumeMaterialSharedPtr();
  }

  // If confined layers exist, loop over them and collect the layer material
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& lay : tVolume.confinedLayers()->arrayObjects()) {
      collectMaterial(*lay, detMatMap);
    }
  }

  // If any of the boundary surfaces has material collect that
  if (m_cfg.processBoundaries) {
    for (auto& bou : tVolume.boundarySurfaces()) {
      const auto& bSurface = bou->surfaceRepresentation();
      if (bSurface.surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[bSurface.geoID()] = bSurface.surfaceMaterialSharedPtr();
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

void FW::RootMaterialWriter::collectMaterial(
    const Acts::Layer& tLayer, Acts::DetectorMaterialMaps& detMatMap) {
  // If the representing surface has material, collect it
  const auto& rSurface = tLayer.surfaceRepresentation();
  if (rSurface.surfaceMaterialSharedPtr() != nullptr and
      m_cfg.processRepresenting) {
    detMatMap.first[rSurface.geoID()] = rSurface.surfaceMaterialSharedPtr();
  }

  // Check the approach surfaces
  if (tLayer.approachDescriptor() != nullptr and m_cfg.processApproaches) {
    for (auto& aSurface : tLayer.approachDescriptor()->containedSurfaces()) {
      if (aSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[aSurface->geoID()] =
            aSurface->surfaceMaterialSharedPtr();
      }
    }
  }

  // Check the sensitive surfaces
  if (tLayer.surfaceArray() != nullptr and m_cfg.processSensitives) {
    // sensitive surface loop
    for (auto& sSurface : tLayer.surfaceArray()->surfaces()) {
      if (sSurface->surfaceMaterialSharedPtr() != nullptr) {
        detMatMap.first[sSurface->geoID()] =
            sSurface->surfaceMaterialSharedPtr();
      }
    }
  }
}
