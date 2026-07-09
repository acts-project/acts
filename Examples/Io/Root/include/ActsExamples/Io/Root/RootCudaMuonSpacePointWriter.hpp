// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#ifdef ACTS_ENABLE_CUDA

#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @brief Root Writer into CudaMuonSpacePointContainer
/// Analogous to RootMuonSpacePointWriter
class RootCudaMuonSpacePointWriter
    : public WriterT<CudaMuonSpacePointContainer> {
 public:
  struct Config {
    /// Input CUDA muon space point collection.
    std::string inputSpacePoints{};

    /// Path to the output ROOT file.
    std::string filePath{};

    /// Output file access mode.
    std::string fileMode{"RECREATE"};

    /// Name of the tree within the output ROOT file.
    std::string treeName{"muonSpacePoints"};
  };

  /// Construct the CUDA muon space point writer.
  ///
  /// @param config The configuration object.
  /// @param level The logging level.
  RootCudaMuonSpacePointWriter(const Config& config,
                               Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootCudaMuonSpacePointWriter() override;

  /// End-of-run hook.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters.
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param ctx The algorithm context.
  /// @param hits The CUDA muon space point container to write.
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const CudaMuonSpacePointContainer& hits) override;

  Config m_cfg{};

  std::unique_ptr<TFile> m_file{};
  TTree* m_tree{};

  mutable std::mutex m_mutex{};

  /// Event identifier.
  std::uint32_t m_eventId{0};

  /// Geometry identifier of the associated surface.
  std::vector<Acts::GeometryIdentifier::Value> m_geometryId{};

  /// Identifier of the associated bucket.
  std::vector<std::uint16_t> m_bucketId{};

  /// Muon identifier.
  std::vector<std::uint32_t> m_muonId{};

  /// Position of the measurement.
  std::vector<float> m_localPositionX{};
  std::vector<float> m_localPositionY{};
  std::vector<float> m_localPositionZ{};

  /// Direction of the sensor / wire.
  std::vector<float> m_sensorDirectionTheta{};
  std::vector<float> m_sensorDirectionPhi{};

  /// Vector pointing to the next channel in the same measurement plane.
  std::vector<float> m_toNextSensorTheta{};
  std::vector<float> m_toNextSensorPhi{};

  /// Covariance value along the non-bending direction.
  std::vector<float> m_covLoc0{};

  /// Covariance value along the bending direction.
  std::vector<float> m_covLoc1{};

  /// Time covariance value.
  std::vector<float> m_covT{};

  /// Drift radius of the straw measurements.
  std::vector<float> m_driftR{};

  /// Recorded measurement time.
  std::vector<float> m_time{};
};

}  // namespace ActsExamples

#endif
