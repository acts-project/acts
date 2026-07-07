// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <boost/container/small_vector.hpp>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootSeedEventAnaWriter
///
/// Write out the seed reconstructed by the seeding algorithm in
/// comma-separated-value format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory.
/// Files are named using the following schema
///
///     event000000001-seed.csv
///     event000000002-seed.csv
///
/// and each line in the file corresponds to one seed.
class RootSeedEventAnaWriter : public WriterT<TrackParametersContainer> {
 public:
  struct Config {
    /// Input estimated track parameters collection.
    std::string inputTrackParameters;
    /// Input seed collection.
    std::string inputSimSeeds;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// output directory
    std::string outputDir;

    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "seedsEvent";
    /// The writing mode
    /// Path to the output ROOT file.
    std::string ROOTFileName = "seedsEvent.root";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  explicit RootSeedEventAnaWriter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Ensure underlying file is closed.
  ~RootSeedEventAnaWriter() final;

  /// End-of-run hook
  ProcessCode finalize() final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trackParams are parameters to write
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackParametersContainer& trackParams) override;

 private:
  Config m_cfg;  ///< The config class

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<SeedContainer> m_inputSimSeeds{this, "InputSimSeeds"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  /// @brief Struct for brief seed summary info
  ///
  struct SeedInfo {
    std::size_t seedID = 0;
    SimBarcode particleId;
    float seedPt = -1;
    float seedPhi = 0;
    float seedEta = 0;
    float vertexZ = 0;
    float quality = -1;
    boost::container::small_vector<Acts::Vector3, 3> globalPosition;
    std::vector<bool> haveTime;
    std::vector<float> spTime;
    std::vector<float> spVarR;
    std::vector<float> spVarZ;
    std::vector<float> spVarT;
    float truthDistance = -1;
    std::string seedType = "unknown";
    int seedTypeInt = 0;  // 0: unknown, 1: good, 2: duplicate, 3: fake
    ProtoTrack measurementsID;
  };

  /// ROOT file and tree
  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  /// Event identifier.
  std::uint32_t m_eventId = 0;
  /// Hit surface identifier.
  std::vector<std::uint64_t> m_measurementId_1;
  std::vector<std::uint64_t> m_measurementId_2;
  std::vector<std::uint64_t> m_measurementId_3;
  /// Space point surface identifier.
  std::vector<std::uint64_t> m_geometryId_1;
  std::vector<std::uint64_t> m_geometryId_2;
  std::vector<std::uint64_t> m_geometryId_3;
  std::vector<int> m_have_time_1;
  std::vector<int> m_have_time_2;
  std::vector<int> m_have_time_3;
  /// Global space point 4D position components in mm. init to NaN
  std::vector<float> m_x_1;
  std::vector<float> m_x_2;
  std::vector<float> m_x_3;
  std::vector<float> m_y_1;
  std::vector<float> m_y_2;
  std::vector<float> m_y_3;
  std::vector<float> m_z_1;
  std::vector<float> m_z_2;
  std::vector<float> m_z_3;
  std::vector<float> m_t_1;
  std::vector<float> m_t_2;
  std::vector<float> m_t_3;
  // Global space point 4D position uncertainties
  std::vector<float> m_var_r_1;
  std::vector<float> m_var_r_2;
  std::vector<float> m_var_r_3;
  std::vector<float> m_var_z_1;
  std::vector<float> m_var_z_2;
  std::vector<float> m_var_z_3;
  std::vector<float> m_var_t_1;
  std::vector<float> m_var_t_2;
  std::vector<float> m_var_t_3;
  // Seed vertex position
  std::vector<double> m_z_vertex;
  // Seed quality
  std::vector<float> m_seed_quality;
  // Seed type int
  std::vector<int> m_seedTypeInt;
  // particle ID
  std::vector<std::uint64_t> m_particleId;

  // for seed estimated parameters
  std::vector<float> m_est_loc0;
  std::vector<float> m_est_loc1;
  std::vector<float> m_est_phi;
  std::vector<float> m_est_theta;
  std::vector<float> m_est_qop;
  std::vector<float> m_est_time;
  std::vector<float> m_est_eta;

  // for seed truth parameters
  // std::vector<float> m_truth_loc0;
  // std::vector<float> m_truth_loc1;
  // std::vector<float> m_truth_phi;
  // std::vector<float> m_truth_theta;
  // std::vector<float> m_truth_qop;
  // std::vector<float> m_truth_time;
  // std::vector<float> m_truth_eta;

  // for seed residuals
  // std::vector<float> m_res_loc0;
  // std::vector<float> m_res_loc1;
  // std::vector<float> m_res_phi;
  // std::vector<float> m_res_theta;
  // std::vector<float> m_res_qop;
  // std::vector<float> m_res_time;
  // std::vector<float> m_res_eta;

  //
  // for volume, layer, surface ids
  std::vector<int> m_volumeId_1;
  std::vector<int> m_volumeId_2;
  std::vector<int> m_volumeId_3;
  std::vector<int> m_layerId_1;
  std::vector<int> m_layerId_2;
  std::vector<int> m_layerId_3;
  std::vector<int> m_surfaceId_1;
  std::vector<int> m_surfaceId_2;
  std::vector<int> m_surfaceId_3;
};

}  // namespace ActsExamples
