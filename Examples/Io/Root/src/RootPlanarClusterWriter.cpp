// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootPlanarClusterWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Digitization/DigitizationSourceLink.hpp"
#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootPlanarClusterWriter::RootPlanarClusterWriter(
    const ActsExamples::RootPlanarClusterWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputClusters, "RootPlanarClusterWriter", level),
      m_cfg(config) {
  // inputClusters is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // Set the branches
  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("volume_id", &m_volumeID);
  m_outputTree->Branch("layer_id", &m_layerID);
  m_outputTree->Branch("surface_id", &m_surfaceID);
  m_outputTree->Branch("g_x", &m_x);
  m_outputTree->Branch("g_y", &m_y);
  m_outputTree->Branch("g_z", &m_z);
  m_outputTree->Branch("g_t", &m_t);
  m_outputTree->Branch("l_x", &m_lx);
  m_outputTree->Branch("l_y", &m_ly);
  m_outputTree->Branch("cov_l_x", &m_cov_lx);
  m_outputTree->Branch("cov_l_y", &m_cov_ly);
  m_outputTree->Branch("cell_ID_x", &m_cell_IDx);
  m_outputTree->Branch("cell_ID_y", &m_cell_IDy);
  m_outputTree->Branch("cell_l_x", &m_cell_lx);
  m_outputTree->Branch("cell_l_y", &m_cell_ly);
  m_outputTree->Branch("cell_data", &m_cell_data);
  m_outputTree->Branch("truth_g_x", &m_t_gx);
  m_outputTree->Branch("truth_g_y", &m_t_gy);
  m_outputTree->Branch("truth_g_z", &m_t_gz);
  m_outputTree->Branch("truth_g_t", &m_t_gt);
  m_outputTree->Branch("truth_l_x", &m_t_lx);
  m_outputTree->Branch("truth_l_y", &m_t_ly);
  m_outputTree->Branch("truth_barcode", &m_t_barcode);
}

ActsExamples::RootPlanarClusterWriter::~RootPlanarClusterWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootPlanarClusterWriter::endRun() {
  // Write the tree
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote clusters to tree '" << m_cfg.treeName << "' in '"
                                       << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootPlanarClusterWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>&
        clusters) {
  // retrieve simulated hits
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);
  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the planar clusters in this event
  for (auto [moduleGeoId, moduleClusters] : groupByModule(clusters)) {
    const Acts::Surface* surfacePtr =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (surfacePtr == nullptr) {
      ACTS_ERROR("Could not find surface for " << moduleGeoId);
      return ProcessCode::ABORT;
    }
    const Acts::Surface& surface = *surfacePtr;

    // geometry identification is the same for all clusters on this module
    m_volumeID = moduleGeoId.volume();
    m_layerID = moduleGeoId.layer();
    m_surfaceID = moduleGeoId.sensitive();

    for (const auto& entry : moduleClusters) {
      const Acts::PlanarModuleCluster& cluster = entry.second;
      // local cluster information: position, @todo coveraiance
      auto parameters = cluster.parameters();
      Acts::Vector2 local(parameters[Acts::BoundIndices::eBoundLoc0],
                          parameters[Acts::BoundIndices::eBoundLoc1]);

      /// prepare for calculating the
      Acts::Vector3 mom(1, 1, 1);
      // transform local into global position information
      Acts::Vector3 pos = surface.localToGlobal(ctx.geoContext, local, mom);
      m_x = pos.x();
      m_y = pos.y();
      m_z = pos.z();
      m_t = parameters[2] / Acts::UnitConstants::ns;
      m_lx = local.x();
      m_ly = local.y();
      m_cov_lx = 0.;  // @todo fill in
      m_cov_ly = 0.;  // @todo fill in
      // get the cells and run through them
      const auto& cells = cluster.digitizationCells();
      auto detectorElement =
          dynamic_cast<const Acts::IdentifiedDetectorElement*>(
              surface.associatedDetectorElement());
      for (auto& cell : cells) {
        // cell identification
        m_cell_IDx.push_back(cell.channel0);
        m_cell_IDy.push_back(cell.channel1);
        m_cell_data.push_back(cell.data);
        // for more we need the digitization module
        if ((detectorElement != nullptr) &&
            detectorElement->digitizationModule()) {
          auto digitationModule = detectorElement->digitizationModule();
          const Acts::Segmentation& segmentation =
              digitationModule->segmentation();
          // get the cell positions
          auto cellLocalPosition = segmentation.cellPosition(cell);
          m_cell_lx.push_back(cellLocalPosition.x());
          m_cell_ly.push_back(cellLocalPosition.y());
        }
      }
      // write hit-particle truth association
      // each hit can have multiple particles, e.g. in a dense environment
      const auto& sl = cluster.sourceLink().get<Acts::DigitizationSourceLink>();
      for (auto idx : sl.indices()) {
        auto it = simHits.nth(idx);
        if (it == simHits.end()) {
          ACTS_FATAL("Simulation hit with index " << idx << " does not exist");
          return ProcessCode::ABORT;
        }
        const auto& simHit = *it;

        // local position to be calculated
        Acts::Vector2 lPosition{0., 0.};
        auto lpResult = surface.globalToLocal(ctx.geoContext, simHit.position(),
                                              simHit.unitDirection());
        if (not lpResult.ok()) {
          ACTS_FATAL("Global to local transformation did not succeed.");
          return ProcessCode::ABORT;
        }
        lPosition = lpResult.value();
        // fill the variables
        m_t_gx.push_back(simHit.position().x());
        m_t_gy.push_back(simHit.position().y());
        m_t_gz.push_back(simHit.position().z());
        m_t_gt.push_back(simHit.time());
        m_t_lx.push_back(lPosition.x());
        m_t_ly.push_back(lPosition.y());
        m_t_barcode.push_back(simHit.particleId().value());
      }
      // fill the tree
      m_outputTree->Fill();
      // now reset
      m_cell_IDx.clear();
      m_cell_IDy.clear();
      m_cell_lx.clear();
      m_cell_ly.clear();
      m_cell_data.clear();
      m_t_gx.clear();
      m_t_gy.clear();
      m_t_gz.clear();
      m_t_gt.clear();
      m_t_lx.clear();
      m_t_ly.clear();
      m_t_barcode.clear();
    }
  }
  return ActsExamples::ProcessCode::SUCCESS;
}
