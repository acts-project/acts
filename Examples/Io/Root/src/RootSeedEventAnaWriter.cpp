// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSeedEventAnaWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootSeedEventAnaWriter::RootSeedEventAnaWriter(
    const ActsExamples::RootSeedEventAnaWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT<TrackParametersContainer>(config.inputTrackParameters,
                                        "RootSeedEventAnaWriter", level),
      m_cfg(config) {
  if (m_cfg.inputSimSeeds.empty()) {
    throw std::invalid_argument("Missing space points input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-simulated-hits map input collection");
  }
  if (m_cfg.outputDir.empty()) {
    throw std::invalid_argument("Missing output directory");
  }
  if (m_cfg.ROOTFileName.empty()) {
    throw std::invalid_argument("Missing output ROOT file name");
  }

  m_inputSimSeeds.initialize(m_cfg.inputSimSeeds);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  // open root file and create the tree
  std::string ROOTFilePath = m_cfg.outputDir + "/" + m_cfg.ROOTFileName;

  m_outputFile = TFile::Open(ROOTFilePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + ROOTFilePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // setup the branches
  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("measurement_id_1", &m_measurementId_1);
  m_outputTree->Branch("measurement_id_2", &m_measurementId_2);
  m_outputTree->Branch("measurement_id_3", &m_measurementId_3);
  m_outputTree->Branch("geometry_id_1", &m_geometryId_1);
  m_outputTree->Branch("geometry_id_2", &m_geometryId_2);
  m_outputTree->Branch("geometry_id_3", &m_geometryId_3);
  m_outputTree->Branch("have_time_1", &m_have_time_1);
  m_outputTree->Branch("have_time_2", &m_have_time_2);
  m_outputTree->Branch("have_time_3", &m_have_time_3);
  m_outputTree->Branch("x_1", &m_x_1);
  m_outputTree->Branch("x_2", &m_x_2);
  m_outputTree->Branch("x_3", &m_x_3);
  m_outputTree->Branch("y_1", &m_y_1);
  m_outputTree->Branch("y_2", &m_y_2);
  m_outputTree->Branch("y_3", &m_y_3);
  m_outputTree->Branch("z_1", &m_z_1);
  m_outputTree->Branch("z_2", &m_z_2);
  m_outputTree->Branch("z_3", &m_z_3);
  m_outputTree->Branch("t_1", &m_t_1);
  m_outputTree->Branch("t_2", &m_t_2);
  m_outputTree->Branch("t_3", &m_t_3);
  m_outputTree->Branch("var_r_1", &m_var_r_1);
  m_outputTree->Branch("var_r_2", &m_var_r_2);
  m_outputTree->Branch("var_r_3", &m_var_r_3);
  m_outputTree->Branch("var_z_1", &m_var_z_1);
  m_outputTree->Branch("var_z_2", &m_var_z_2);
  m_outputTree->Branch("var_z_3", &m_var_z_3);
  m_outputTree->Branch("var_t_1", &m_var_t_1);
  m_outputTree->Branch("var_t_2", &m_var_t_2);
  m_outputTree->Branch("var_t_3", &m_var_t_3);
  m_outputTree->Branch("z_vertex", &m_z_vertex);
  m_outputTree->Branch("seed_quality", &m_seed_quality);
  m_outputTree->Branch("seed_type", &m_seedTypeInt);
  m_outputTree->Branch("particle_id", &m_particleId);
  //
  m_outputTree->Branch("est_loc0", &m_est_loc0);
  m_outputTree->Branch("est_loc1", &m_est_loc1);
  m_outputTree->Branch("est_phi", &m_est_phi);
  m_outputTree->Branch("est_theta", &m_est_theta);
  m_outputTree->Branch("est_qop", &m_est_qop);
  m_outputTree->Branch("est_time", &m_est_time);
  m_outputTree->Branch("est_eta", &m_est_eta);

  // for seed truth parameters
  // m_outputTree->Branch("truth_loc0", &m_truth_loc0);
  // m_outputTree->Branch("truth_loc1", &m_truth_loc1);
  // m_outputTree->Branch("truth_phi", &m_truth_phi);
  // m_outputTree->Branch("truth_theta", &m_truth_theta);
  // m_outputTree->Branch("truth_qop", &m_truth_qop);
  // m_outputTree->Branch("truth_time", &m_truth_time);
  // m_outputTree->Branch("truth_eta", &m_truth_eta);
  // for residuals
  // m_outputTree->Branch("res_loc0", &m_res_loc0);
  // m_outputTree->Branch("res_loc1", &m_res_loc1);
  // m_outputTree->Branch("res_phi", &m_res_phi);
  // m_outputTree->Branch("res_theta", &m_res_theta);
  // m_outputTree->Branch("res_qop", &m_res_qop);
  // m_outputTree->Branch("res_time", &m_res_time);
  // m_outputTree->Branch("res_eta", &m_res_eta);
}

ActsExamples::RootSeedEventAnaWriter::~RootSeedEventAnaWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootSeedEventAnaWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_VERBOSE("Wrote seeds to tree '"
               << m_cfg.treeName << "' in '"
               << m_cfg.outputDir + "/" + m_cfg.ROOTFileName << "'");
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootSeedEventAnaWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const TrackParametersContainer& trackParams) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // reset all variables
  m_eventId = 0;
  m_measurementId_1.clear();
  m_measurementId_2.clear();
  m_measurementId_3.clear();
  m_geometryId_1.clear();
  m_geometryId_2.clear();
  m_geometryId_3.clear();
  m_have_time_1.clear();
  m_have_time_2.clear();
  m_have_time_3.clear();
  m_x_1.clear();
  m_x_2.clear();
  m_x_3.clear();
  m_y_1.clear();
  m_y_2.clear();
  m_y_3.clear();
  m_z_1.clear();
  m_z_2.clear();
  m_z_3.clear();
  m_t_1.clear();
  m_t_2.clear();
  m_t_3.clear();
  m_var_r_1.clear();
  m_var_r_2.clear();
  m_var_r_3.clear();
  m_var_z_1.clear();
  m_var_z_2.clear();
  m_var_z_3.clear();
  m_var_t_1.clear();
  m_var_t_2.clear();
  m_var_t_3.clear();
  m_z_vertex.clear();
  m_seed_quality.clear();
  m_seedTypeInt.clear();
  m_particleId.clear();

  m_est_loc0.clear();
  m_est_loc1.clear();
  m_est_phi.clear();
  m_est_theta.clear();
  m_est_qop.clear();
  m_est_time.clear();
  m_est_eta.clear();
  // m_truth_loc0.clear();
  // m_truth_loc1.clear();
  // m_truth_phi.clear();
  // m_truth_theta.clear();
  // m_truth_qop.clear();
  // m_truth_time.clear();
  // m_truth_eta.clear();
  // m_res_loc0.clear();
  // m_res_loc1.clear();
  // m_res_phi.clear();
  // m_res_theta.clear();
  // m_res_qop.clear();
  // m_res_time.clear();
  // m_res_eta.clear();

  // Read additional input collections
  const auto& seeds = m_inputSimSeeds(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // Get the event number
  m_eventId = ctx.eventNumber;

  std::unordered_map<std::size_t, SeedInfo> infoMap;
  std::unordered_map<SimBarcode, std::pair<std::size_t, float>> goodSeed;

  // Loop over the estimated track parameters
  for (std::size_t iparams = 0; iparams < trackParams.size(); ++iparams) {
    // The estimated bound parameters vector
    const auto& seedEstParams = trackParams[iparams];
    const auto params = trackParams[iparams].parameters();

    // The reference surface of the parameters, i.e. also the reference surface
    // of the first space point
    // const auto& surface = params.referenceSurface();
    // m_volumeId = surface.geometryId().volume();
    // m_layerId = surface.geometryId().layer();
    // m_surfaceId = surface.geometryId().sensitive();

    // Get and fill the track estimated parameters
    float seedPhi = params[Acts::eBoundPhi];
    float seedEta = std::atanh(std::cos(params[Acts::eBoundTheta]));
    m_est_loc0.push_back(params[Acts::eBoundLoc0]);
    m_est_loc1.push_back(params[Acts::eBoundLoc1]);
    m_est_phi.push_back(params[Acts::eBoundPhi]);
    m_est_theta.push_back(params[Acts::eBoundTheta]);
    m_est_qop.push_back(params[Acts::eBoundQOverP]);
    m_est_time.push_back(params[Acts::eBoundTime]);

    ACTS_VERBOSE("eta: " << eta(seedEstParams.momentum())
                         << ", seedEta: " << seedEta);

    // Get the proto track from which the track parameters are estimated
    const auto& seed = seeds[iparams];
    const auto& ptrack = seedToProtoTrack(seed);

    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, ptrack, particleHitCounts);
    bool truthMatched = false;
    float truthDistance = -1;
    auto majorityParticleId = particleHitCounts.front().particleId;
    // Seed are considered truth matched if they have only one contributing
    // particle
    if (particleHitCounts.size() == 1) {
      truthMatched = true;
      // Get the index of the first space point
      const auto& hitIdx = ptrack.front();
      // Get the sim hits via the measurement to sim hits map
      auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
      // Get the truth particle direction from the sim hits
      Acts::Vector3 truthUnitDir = {0, 0, 0};
      for (auto [_, simHitIdx] : indices) {
        const auto& simHit = *simHits.nth(simHitIdx);
        if (simHit.particleId() == majorityParticleId) {
          truthUnitDir = simHit.direction();
        }
      }
      // Compute the distance between the truth and estimated directions
      float truthPhi = phi(truthUnitDir);
      float truthEta = std::atanh(std::cos(theta(truthUnitDir)));
      float dEta = std::abs(truthEta - seedEta);
      float dPhi =
          std::abs(truthPhi - seedPhi) < std::numbers::pi_v<float>
              ? std::abs(truthPhi - seedPhi)
              : std::abs(truthPhi - seedPhi) - std::numbers::pi_v<float>;
      truthDistance = std::sqrt(dPhi * dPhi + dEta * dEta);
      // If the seed is truth matched, check if it is the closest one for the
      // contributing particle
      if (goodSeed.contains(majorityParticleId)) {
        if (goodSeed[majorityParticleId].second > truthDistance) {
          goodSeed[majorityParticleId] = std::make_pair(iparams, truthDistance);
        }
      } else {
        goodSeed[majorityParticleId] = std::make_pair(iparams, truthDistance);
      }
    }
    // Store the global position of the space points
    boost::container::small_vector<Acts::Vector3, 3> globalPosition;
    std::vector<bool> haveTime;
    std::vector<float> spTime;
    std::vector<float> spVarR;
    std::vector<float> spVarZ;
    std::vector<float> spVarT;
    for (auto spacePointPtr : seed.spacePoints()) {
      Acts::Vector3 pos(spacePointPtr.x(), spacePointPtr.y(),
                        spacePointPtr.z());
      globalPosition.push_back(pos);
      if (std::isnan(spacePointPtr.time())) {
        haveTime.push_back(false);
      } else {
        haveTime.push_back(true);
      }
      spTime.push_back(spacePointPtr.time());
      spVarR.push_back(spacePointPtr.varianceR());
      spVarZ.push_back(spacePointPtr.varianceZ());
      spVarT.push_back(spacePointPtr.varianceT());
    }

    // track info
    SeedInfo toAdd;
    toAdd.seedID = iparams;
    toAdd.particleId = majorityParticleId;
    toAdd.seedPt = std::abs(1.0 / params[Acts::eBoundQOverP]) *
                   std::sin(params[Acts::eBoundTheta]);
    toAdd.seedPhi = seedPhi;
    toAdd.seedEta = seedEta;
    toAdd.vertexZ = seed.vertexZ();
    toAdd.quality = seed.quality();
    toAdd.globalPosition = globalPosition;
    toAdd.haveTime = haveTime;
    toAdd.spTime = spTime;
    toAdd.spVarR = spVarR;
    toAdd.spVarZ = spVarZ;
    toAdd.spVarT = spVarT;
    toAdd.truthDistance = truthDistance;
    toAdd.seedType = truthMatched ? "duplicate" : "fake";
    toAdd.seedTypeInt = truthMatched ? 2 : 3;
    toAdd.measurementsID = ptrack;

    infoMap[toAdd.seedID] = toAdd;
  }

  for (auto& [id, info] : infoMap) {
    if (goodSeed[info.particleId].first == id) {
      info.seedType = "good";
      info.seedTypeInt = 1;
    }

    // Fill the tree
    m_measurementId_1.push_back(info.measurementsID[0]);
    m_measurementId_2.push_back(info.measurementsID[1]);
    m_measurementId_3.push_back(info.measurementsID[2]);

    m_x_1.push_back(info.globalPosition[0].x());
    m_x_2.push_back(info.globalPosition[1].x());
    m_x_3.push_back(info.globalPosition[2].x());
    m_y_1.push_back(info.globalPosition[0].y());
    m_y_2.push_back(info.globalPosition[1].y());
    m_y_3.push_back(info.globalPosition[2].y());
    m_z_1.push_back(info.globalPosition[0].z());
    m_z_2.push_back(info.globalPosition[1].z());
    m_z_3.push_back(info.globalPosition[2].z());
    m_have_time_1.push_back(info.haveTime[0]);
    m_have_time_2.push_back(info.haveTime[1]);
    m_have_time_3.push_back(info.haveTime[2]);
    m_t_1.push_back(info.spTime[0]);
    m_t_2.push_back(info.spTime[1]);
    m_t_3.push_back(info.spTime[2]);
    m_var_r_1.push_back(info.spVarR[0]);
    m_var_r_2.push_back(info.spVarR[1]);
    m_var_r_3.push_back(info.spVarR[2]);
    m_var_z_1.push_back(info.spVarZ[0]);
    m_var_z_2.push_back(info.spVarZ[1]);
    m_var_z_3.push_back(info.spVarZ[2]);
    m_var_t_1.push_back(info.spVarT[0]);
    m_var_t_2.push_back(info.spVarT[1]);
    m_var_t_3.push_back(info.spVarT[2]);
    m_z_vertex.push_back(info.vertexZ);
    m_seed_quality.push_back(info.quality);
    m_seedTypeInt.push_back(info.seedTypeInt);
    m_particleId.push_back(info.particleId.hash());
  }
  m_outputTree->Fill();
  return ProcessCode::SUCCESS;
}
