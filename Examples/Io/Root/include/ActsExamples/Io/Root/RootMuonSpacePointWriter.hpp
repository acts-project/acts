// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

class RootMuonSpacePointWriter : public WriterT<MuonSpacePointContainer> {
  struct Config {
    /// Input sim hit collection to write.
    std::string inputSpacePoints{};
    /// Path to the output file.
    std::string filePath{};
    /// Output file access mode.
    std::string fileMode{"RECREATE"};
    /// Name of the tree within the output file.
    std::string treeName{"muonSpacePoints"};
  };

  /// Construct the particle writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootMuonSpacePointWriter(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootMuonSpacePointWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] hits are the hits to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MuonSpacePointContainer& hits) override;

  Config m_cfg{};
  std::unique_ptr<TFile> m_file{};
  TTree* m_tree{};

  /// @brief Event identifier.
  std::uint32_t m_eventId{0};
  /// @brief Geometry identifier of the associated surface
  std::vector<Acts::GeometryIdentifier::Value> m_geometryId{};
  /// @brief Position of the measurement
  std::vector<float> m_localPositionX{};
  std::vector<float> m_localPositionY{};
  std::vector<float> m_localPositionZ{};
  /// @brief Direction of the sensor / wire
  std::vector<float> m_sensorDirectionX{};
  std::vector<float> m_sensorDirectionY{};
  std::vector<float> m_sensorDirectionZ{};
  /// @brief Vector pointing to the next channel in the same measurement plane
  std::vector<float> m_toNextSensorX{};
  std::vector<float> m_toNextSensorY{};
  std::vector<float> m_toNextSensorZ{};
  /// @brief Covariance value along the non-bending direction
  std::vector<float> m_covLoc0{};
  /// @brief Covaraiance value along the bending direction
  std::vector<float> m_covLoc1{};
  /// @brief Time covariance value
  std::vector<float> m_covLocT{};
  /// @brief Drift radius of the straw measurements
  std::vector<float> m_driftR{};
  /// @brief Recorded measurement time.
  std::vector<float> m_time{};
};

}  // namespace ActsExamples
