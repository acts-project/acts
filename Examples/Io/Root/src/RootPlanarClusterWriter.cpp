// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootPlanarClusterWriter.hpp"

#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>

#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimIdentifier.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Digitization/Segmentation.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Utilities/Units.hpp"

FW::RootPlanarClusterWriter::RootPlanarClusterWriter(
    const FW::RootPlanarClusterWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputClusters, "RootPlanarClusterWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // inputClusters is already checked by base constructor
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.treeName.c_str(), "TTree from RootPlanarClusterWriter");
  if (m_outputTree == nullptr)
    throw std::bad_alloc();

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
  m_outputTree->Branch("truth_barcode", &m_t_barcode, "truth_barcode/l");
}

FW::RootPlanarClusterWriter::~RootPlanarClusterWriter() {
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::RootPlanarClusterWriter::endRun() {
  // Write the tree
  m_outputFile->cd();
  m_outputTree->Write();
  ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                        << m_cfg.filePath << "'");
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::RootPlanarClusterWriter::writeT(
    const AlgorithmContext& ctx,
    const FW::GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters) {
  // retrieve simulated hits
  const auto& simHits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);
  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the planar clusters in this event
  for (const auto& entry : clusters) {
    Acts::GeometryID geoId = entry.first;
    const Acts::PlanarModuleCluster& cluster = entry.second;
    // local cluster information: position, @todo coveraiance
    auto parameters = cluster.parameters();
    Acts::Vector2D local(parameters[Acts::ParDef::eLOC_0],
                         parameters[Acts::ParDef::eLOC_1]);

    /// prepare for calculating the
    Acts::Vector3D pos(0, 0, 0);
    Acts::Vector3D mom(1, 1, 1);
    // the cluster surface
    const auto& clusterSurface = cluster.referenceSurface();
    // transform local into global position information
    clusterSurface.localToGlobal(ctx.geoContext, local, mom, pos);
    // identification
    m_volumeID = geoId.volume();
    m_layerID = geoId.layer();
    m_surfaceID = geoId.sensitive();
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
    auto detectorElement = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
        clusterSurface.associatedDetectorElement());
    for (auto& cell : cells) {
      // cell identification
      m_cell_IDx.push_back(cell.channel0);
      m_cell_IDy.push_back(cell.channel1);
      m_cell_data.push_back(cell.data);
      // for more we need the digitization module
      if (detectorElement && detectorElement->digitizationModule()) {
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
    for (auto idx : cluster.sourceLink().indices()) {
      auto it = simHits.nth(idx);
      if (it == simHits.end()) {
        ACTS_FATAL("Simulation hit with index " << idx << " does not exist");
        return ProcessCode::ABORT;
      }
      const auto& simHit = *it;

      // local position to be calculated
      Acts::Vector2D lPosition;
      clusterSurface.globalToLocal(ctx.geoContext, simHit.position(),
                                   simHit.unitDirection(), lPosition);
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
  return FW::ProcessCode::SUCCESS;
}
