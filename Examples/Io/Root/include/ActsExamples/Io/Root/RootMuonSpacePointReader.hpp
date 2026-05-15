// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"

namespace ActsExamples {

class RootMuonSpacePointReader : public IReader {
 public:
  struct Config {
    /// Input sim hit collection to write.
    std::string outputSpacePoints{"MuonSpacePoints"};
    /// Path to the output file.
    std::string filePath{};
    /// Name of the tree within the output file.
    std::string treeName{"muonSpacePoints"};
  };

  /// Construct the particle writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootMuonSpacePointReader(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootMuonSpacePointReader() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  /// Framework name() method
  std::string name() const override { return logger().name(); }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const AlgorithmContext& context) override;

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

  mutable std::mutex m_mutex{};

  /// @brief Configuration object
  Config m_cfg{};

  WriteDataHandle<MuonSpacePointContainer> m_outputContainer{
      this, "OutputSpacePoints"};
  std::unique_ptr<const Acts::Logger> m_logger{};

  /// @brief Input file directly read at construction stage
  std::unique_ptr<TFile> m_file{TFile::Open(m_cfg.filePath.c_str(), "READ")};

  /// @brief TTree reader
  TTreeReader m_reader{m_cfg.treeName.c_str(), m_file.get()};

  std::vector<std::uint32_t> m_eventRanges{};

  template <typename T>
  using VecReader_t = TTreeReaderValue<std::vector<T>>;
  /// @brief Event identifier.
  TTreeReaderValue<std::uint32_t> m_eventId{m_reader, "event_id"};

  /// @brief Geometry identifier of the associated surface
  VecReader_t<Acts::GeometryIdentifier::Value> m_geometryId{
      m_reader, "spacePoint_geometryId"};
  /// @brief Identifier of the associated bucket
  VecReader_t<std::uint16_t> m_bucketId{m_reader, "spacePoint_bucketId"};
  /// @brief Muon identifier
  VecReader_t<std::uint32_t> m_muonId{m_reader, "spacePoint_muonId"};
  /// @brief Position of the measurement
  VecReader_t<float> m_localPositionX{m_reader, "spacePoint_localPosX"};
  VecReader_t<float> m_localPositionY{m_reader, "spacePoint_localPosY"};
  VecReader_t<float> m_localPositionZ{m_reader, "spacePoint_localPosZ"};
  /// @brief Direction of the sensor / wire
  VecReader_t<float> m_sensorDirectionTheta{m_reader,
                                            "spacePoint_sensorDirTheta"};
  VecReader_t<float> m_sensorDirectionPhi{m_reader, "spacePoint_sensorDirPhi"};
  /// @brief Vector pointing to the next channel in the same measurement plane
  VecReader_t<float> m_toNextSensorTheta{m_reader, "spacePoint_toNextDirTheta"};
  VecReader_t<float> m_toNextSensorPhi{m_reader, "spacePoint_toNextDirPhi"};
  /// @brief Covariance value along the non-bending direction
  VecReader_t<float> m_covLoc0{m_reader, "spacePoint_covLoc0"};
  /// @brief Covaraiance value along the bending direction
  VecReader_t<float> m_covLoc1{m_reader, "spacePoint_covLoc1"};
  /// @brief Time covariance value
  VecReader_t<float> m_covT{m_reader, "spacePoint_covT"};
  /// @brief Drift radius of the straw measurements
  VecReader_t<float> m_driftR{m_reader, "spacePoint_driftRadius"};
  /// @brief Recorded measurement time.
  VecReader_t<float> m_time{m_reader, "spacePoint_time"};
};

}  // namespace ActsExamples
