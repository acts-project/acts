// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ios>
#include <stdexcept>
#include <unordered_map>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootAthenaDumpWriter::RootAthenaDumpWriter(const Config& config,
                                           Acts::Logging::Level level)
    : IWriter(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputClusters.initialize(m_cfg.inputClusters);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputMeasParticleMap.initialize(m_cfg.inputMeasParticleMap);
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);

  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), "RECREATE");
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // Reserve capacity upfront so that data() pointers passed to ROOT branches
  // remain stable for the lifetime of the writer. The vectors are cleared and
  // refilled each event; clear() preserves capacity, keeping the pointer valid.
  m_partEventNumber.reserve(s_maxParticles);
  m_partBarcode.reserve(s_maxParticles);
  m_partPt.reserve(s_maxParticles);
  m_partEta.reserve(s_maxParticles);
  m_partPdgId.reserve(s_maxParticles);
  m_partVx.reserve(s_maxParticles);
  m_partVy.reserve(s_maxParticles);
  m_partVz.reserve(s_maxParticles);

  m_clIndex.reserve(s_maxClusters);
  m_clBarrelEndcap.reserve(s_maxClusters);
  m_clEtaModule.reserve(s_maxClusters);
  m_clPhiModule.reserve(s_maxClusters);
  m_clModuleId.reserve(s_maxClusters);
  m_clX.reserve(s_maxClusters);
  m_clY.reserve(s_maxClusters);
  m_clZ.reserve(s_maxClusters);

  m_spIndex.reserve(s_maxSpacePoints);
  m_spX.reserve(s_maxSpacePoints);
  m_spY.reserve(s_maxSpacePoints);
  m_spZ.reserve(s_maxSpacePoints);
  m_spCL1Index.reserve(s_maxSpacePoints);
  m_spCL2Index.reserve(s_maxSpacePoints);
  m_spIsOverlap.reserve(s_maxSpacePoints);

  // Event scalar
  m_outputTree->Branch("event_number", &m_eventNumber);

  // Particle branches – variable-length C-arrays, count driven by nPartEVT
  m_outputTree->Branch("nPartEVT", &m_nPartEVT, "nPartEVT/I");
  m_outputTree->Branch("Part_event_number", m_partEventNumber.data(),
                       "Part_event_number[nPartEVT]/I");
  m_outputTree->Branch("Part_barcode", m_partBarcode.data(),
                       "Part_barcode[nPartEVT]/I");
  m_outputTree->Branch("Part_pt", m_partPt.data(), "Part_pt[nPartEVT]/F");
  m_outputTree->Branch("Part_eta", m_partEta.data(), "Part_eta[nPartEVT]/F");
  m_outputTree->Branch("Part_pdg_id", m_partPdgId.data(),
                       "Part_pdg_id[nPartEVT]/I");
  m_outputTree->Branch("Part_vx", m_partVx.data(), "Part_vx[nPartEVT]/F");
  m_outputTree->Branch("Part_vy", m_partVy.data(), "Part_vy[nPartEVT]/F");
  m_outputTree->Branch("Part_vz", m_partVz.data(), "Part_vz[nPartEVT]/F");

  // Cluster branches
  m_outputTree->Branch("nCL", &m_nCL, "nCL/I");
  m_outputTree->Branch("CLindex", m_clIndex.data(), "CLindex[nCL]/I");
  m_outputTree->Branch("CLhardware", &m_clHardware);
  m_outputTree->Branch("CLbarrel_endcap", m_clBarrelEndcap.data(),
                       "CLbarrel_endcap[nCL]/I");
  m_outputTree->Branch("CLeta_module", m_clEtaModule.data(),
                       "CLeta_module[nCL]/I");
  m_outputTree->Branch("CLphi_module", m_clPhiModule.data(),
                       "CLphi_module[nCL]/I");
  m_outputTree->Branch("CLmoduleID", m_clModuleId.data(), "CLmoduleID[nCL]/l");
  m_outputTree->Branch("CLx", m_clX.data(), "CLx[nCL]/D");
  m_outputTree->Branch("CLy", m_clY.data(), "CLy[nCL]/D");
  m_outputTree->Branch("CLz", m_clZ.data(), "CLz[nCL]/D");
  m_outputTree->Branch("CLlocal_cov", &m_clLocalCov);
  m_outputTree->Branch("CLparticleLink_eventIndex",
                       &m_clParticleLinkEventIndex);
  m_outputTree->Branch("CLparticleLink_barcode", &m_clParticleLinkBarcode);

  // Space point branches
  m_outputTree->Branch("nSP", &m_nSP, "nSP/I");
  m_outputTree->Branch("SPindex", m_spIndex.data(), "SPindex[nSP]/I");
  m_outputTree->Branch("SPx", m_spX.data(), "SPx[nSP]/D");
  m_outputTree->Branch("SPy", m_spY.data(), "SPy[nSP]/D");
  m_outputTree->Branch("SPz", m_spZ.data(), "SPz[nSP]/D");
  m_outputTree->Branch("SPCL1_index", m_spCL1Index.data(),
                       "SPCL1_index[nSP]/I");
  m_outputTree->Branch("SPCL2_index", m_spCL2Index.data(),
                       "SPCL2_index[nSP]/I");
  m_outputTree->Branch("SPisOverlap", m_spIsOverlap.data(),
                       "SPisOverlap[nSP]/I");
}

RootAthenaDumpWriter::~RootAthenaDumpWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootAthenaDumpWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();
  ACTS_VERBOSE("Wrote Athena dump to '" << m_cfg.filePath << "'");
  return ProcessCode::SUCCESS;
}

ProcessCode RootAthenaDumpWriter::write(const AlgorithmContext& ctx) {
  const auto& particles = m_inputParticles(ctx);
  const auto& clusters = m_inputClusters(ctx);
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& measPartMap = m_inputMeasParticleMap(ctx);
  const auto& spacePoints = m_inputSpacePoints(ctx);

  std::lock_guard<std::mutex> lock(m_writeMutex);

  if (particles.size() > s_maxParticles) {
    throw std::runtime_error(
        "RootAthenaDumpWriter: particle count exceeds maxParticles (" +
        std::to_string(particles.size()) + " > " +
        std::to_string(s_maxParticles) + ")");
  }
  if (measurements.size() > s_maxClusters) {
    throw std::runtime_error(
        "RootAthenaDumpWriter: measurement count exceeds maxClusters (" +
        std::to_string(measurements.size()) + " > " +
        std::to_string(s_maxClusters) + ")");
  }
  if (spacePoints.size() > s_maxSpacePoints) {
    throw std::runtime_error(
        "RootAthenaDumpWriter: space point count exceeds maxSpacePoints (" +
        std::to_string(spacePoints.size()) + " > " +
        std::to_string(s_maxSpacePoints) + ")");
  }

  m_eventNumber = ctx.eventNumber;

  // Build barcode map: assign sequential Athena barcodes per event.
  // Primaries (vertexSecondary == 0) get barcodes [1, s_maxBarcodeForPrimary].
  // Secondaries get barcodes starting at s_maxBarcodeForPrimary + 1.
  // The same values are used in the particle table and cluster particle links
  // so that make_particle_id(subevent, barcode) is consistent throughout.
  std::unordered_map<ActsFatras::Barcode, int> barcodeMap;
  barcodeMap.reserve(particles.size());
  int primaryCount = 1;
  int secondaryCount = s_maxBarcodeForPrimary + 1;
  for (const auto& particle : particles) {
    const bool isPrimary = particle.particleId().vertexSecondary() == 0;
    barcodeMap[particle.particleId()] =
        isPrimary ? primaryCount++ : secondaryCount++;
  }

  // --- Particles ---
  m_partEventNumber.clear();
  m_partBarcode.clear();
  m_partPt.clear();
  m_partEta.clear();
  m_partPdgId.clear();
  m_partVx.clear();
  m_partVy.clear();
  m_partVz.clear();

  for (const auto& particle : particles) {
    const float pT = static_cast<float>(particle.transverseMomentum() /
                                        Acts::UnitConstants::MeV);
    const float eta = static_cast<float>(
        Acts::VectorHelpers::eta(particle.momentum().normalized()));

    m_partEventNumber.push_back(
        static_cast<int>(particle.particleId().vertexPrimary()));
    m_partBarcode.push_back(barcodeMap.at(particle.particleId()));
    m_partPt.push_back(pT);
    m_partEta.push_back(eta);
    m_partPdgId.push_back(static_cast<int>(particle.pdg()));
    m_partVx.push_back(
        static_cast<float>(particle.position().x() / Acts::UnitConstants::mm));
    m_partVy.push_back(
        static_cast<float>(particle.position().y() / Acts::UnitConstants::mm));
    m_partVz.push_back(
        static_cast<float>(particle.position().z() / Acts::UnitConstants::mm));
  }
  m_nPartEVT = static_cast<int>(m_partBarcode.size());

  // --- Clusters ---
  m_clIndex.clear();
  m_clHardware.clear();
  m_clBarrelEndcap.clear();
  m_clEtaModule.clear();
  m_clPhiModule.clear();
  m_clModuleId.clear();
  m_clX.clear();
  m_clY.clear();
  m_clZ.clear();
  m_clLocalCov.clear();
  m_clParticleLinkEventIndex.clear();
  m_clParticleLinkBarcode.clear();

  for (std::size_t im = 0; im < measurements.size(); ++im) {
    const auto meas = measurements.at(im);

    const bool isPixel = (meas.size() == 2);

    m_clIndex.push_back(static_cast<int>(im));
    m_clHardware.push_back(isPixel ? "PIXEL" : "STRIP");
    m_clModuleId.push_back(meas.geometryId().value());

    m_clBarrelEndcap.push_back(s_sentinel);
    m_clEtaModule.push_back(s_sentinel);
    m_clPhiModule.push_back(s_sentinel);

    const Acts::Vector3& pos = clusters.at(im).globalPosition;
    m_clX.push_back(pos.x());
    m_clY.push_back(pos.y());
    m_clZ.push_back(pos.z());

    if (isPixel) {
      m_clLocalCov.push_back({meas.covariance()(0, 0), meas.covariance()(0, 1),
                              meas.covariance()(1, 0),
                              meas.covariance()(1, 1)});
    } else {
      const double var = meas.covariance()(0, 0);
      m_clLocalCov.push_back({var, var});
    }

    std::vector<int> eventIndices;
    std::vector<int> barcodes;
    const auto [begin, end] = measPartMap.equal_range(im);
    for (auto it = begin; it != end; ++it) {
      const auto& bc = it->second;
      eventIndices.push_back(s_sentinel);
      const auto bcIt = barcodeMap.find(bc);
      barcodes.push_back(bcIt != barcodeMap.end() ? bcIt->second : 0);
    }
    m_clParticleLinkEventIndex.push_back(std::move(eventIndices));
    m_clParticleLinkBarcode.push_back(std::move(barcodes));
  }
  m_nCL = static_cast<int>(m_clIndex.size());

  // --- Space points ---
  m_spIndex.clear();
  m_spX.clear();
  m_spY.clear();
  m_spZ.clear();
  m_spCL1Index.clear();
  m_spCL2Index.clear();
  m_spIsOverlap.clear();

  for (const auto& sp : spacePoints) {
    const auto sLinks = sp.sourceLinks();
    if (sLinks.empty()) {
      ACTS_WARNING("Space point with no source links, skipping");
      continue;
    }

    m_spIndex.push_back(static_cast<int>(m_spIndex.size()));
    m_spX.push_back(static_cast<double>(sp.x()));
    m_spY.push_back(static_cast<double>(sp.y()));
    m_spZ.push_back(static_cast<double>(sp.z()));

    m_spCL1Index.push_back(
        static_cast<int>(sLinks[0].get<IndexSourceLink>().index()));
    m_spCL2Index.push_back(
        sLinks.size() > 1
            ? static_cast<int>(sLinks[1].get<IndexSourceLink>().index())
            : -1);

    // ACTS simulation does not produce overlapping SPs
    m_spIsOverlap.push_back(0);
  }
  m_nSP = static_cast<int>(m_spIndex.size());

  m_outputTree->Fill();

  ACTS_DEBUG("Wrote event " << m_eventNumber << " with " << m_nPartEVT
                            << " particles, " << m_nCL << " clusters, and "
                            << m_nSP << " space points");

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
