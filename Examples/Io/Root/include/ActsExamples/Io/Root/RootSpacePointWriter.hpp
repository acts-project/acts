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
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <limits>
#include <mutex>
#include <string>

class TFile;
class TTree;

namespace ActsExamples {

/// Write out space points as a flat TTree.
///
/// Each entry in the TTree corresponds to one space point for optimum writing
/// speed. The event number is part of the written data.
///
/// Safe to use from multiple writer threads. To avoid thread-saftey issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootSpacePointWriter final : public WriterT<SpacePointContainer> {
 public:
  struct Config {
    /// Input particle collection to write.
    std::string inputSpacePoints;
    /// Input sim hits collection (optional, only needed for residuals)
    std::string inputSimHits;
    /// Input measurement particles map (optional)
    std::string inputMeasurementParticlesMap;
    /// Input measurement to sim hits map (optional, only needed for residuals)
    std::string inputMeasurementSimHitsMap;

    /// Tracking geometry for transformation lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "spacePoints";
  };

  /// Construct the particle writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootSpacePointWriter(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootSpacePointWriter() final;

  /// End-of-run hook
  ProcessCode finalize() final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] spacePoints are the space points to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SpacePointContainer& spacePoints) final;

 private:
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  Config m_cfg;

  std::mutex m_writeMutex;

  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;

  /// Event identifier.
  std::uint32_t m_eventId = 0;
  /// Hit surface identifier.
  std::uint64_t m_measurementId1 = 0;
  std::uint64_t m_measurementId2 = 0;
  /// Space point surface identifier.
  std::uint64_t m_geometryId1 = 0;
  std::uint32_t m_volumeId1 = 0;
  std::uint32_t m_layerId1 = 0;
  std::uint32_t m_surfaceId1 = 0;
  std::uint32_t m_extraId1 = 0;
  std::uint64_t m_geometryId2 = 0;
  std::uint32_t m_volumeId2 = 0;
  std::uint32_t m_layerId2 = 0;
  std::uint32_t m_surfaceId2 = 0;
  std::uint32_t m_extraId2 = 0;
  /// Global space point position components in mm.
  float m_x = std::numeric_limits<float>::quiet_NaN();
  float m_y = std::numeric_limits<float>::quiet_NaN();
  float m_z = std::numeric_limits<float>::quiet_NaN();
  float m_t = std::numeric_limits<float>::quiet_NaN();
  float m_r = std::numeric_limits<float>::quiet_NaN();
  // Global space point position uncertainties
  float m_var_r = std::numeric_limits<float>::quiet_NaN();
  float m_var_z = std::numeric_limits<float>::quiet_NaN();

  // Fake space point (only relevant for strip)
  bool m_fake{};

  float m_trueX = std::numeric_limits<float>::quiet_NaN();
  float m_trueY = std::numeric_limits<float>::quiet_NaN();
  float m_trueZ = std::numeric_limits<float>::quiet_NaN();
  float m_trueT = std::numeric_limits<float>::quiet_NaN();
  float m_trueR = std::numeric_limits<float>::quiet_NaN();

  float m_residualRPhi = std::numeric_limits<float>::quiet_NaN();
  float m_residualZ = std::numeric_limits<float>::quiet_NaN();
  float m_residualDCA = std::numeric_limits<float>::quiet_NaN();
};

}  // namespace ActsExamples
