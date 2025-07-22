// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cstddef>
#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootMaterialTrackWriter::RootMaterialTrackWriter(
    const RootMaterialTrackWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMaterialTracks, "RootMaterialTrackWriter", level),
      m_cfg(config),
      m_accessor({config.prePostStep, config.storeSurface, config.storeVolume,
                  config.recalculateTotals}) {
  // An input collection name and tree name must be specified
  if (m_cfg.inputMaterialTracks.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.treeName.c_str(), "TTree from RootMaterialTrackWriter");
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }
  // Connect the branches
  m_accessor.connectForWrite(*m_outputTree);
}

RootMaterialTrackWriter::~RootMaterialTrackWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootMaterialTrackWriter::finalize() {
  // write the tree and close the file
  ACTS_INFO("Writing ROOT output File : " << m_cfg.filePath);

  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootMaterialTrackWriter::writeT(
    const AlgorithmContext& ctx,
    const std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>&
        materialTracks) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);
  // Loop over the material tracks and write them out
  for (auto& [idTrack, mtrack] : materialTracks) {
    // write & fill
    m_accessor.write(ctx.geoContext, ctx.eventNumber, mtrack);
    m_outputTree->Fill();
  }
  // return success
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
