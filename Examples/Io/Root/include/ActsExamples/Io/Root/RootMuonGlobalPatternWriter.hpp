// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/TrackFinding/GlobalPatternFinderAlgorithm.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @brief Writes the patterns found by the GlobalPatternFinderAlgorithm into a
///        n-tuple. Each entry holds all patterns of one event: the pattern
///        parameters and, for every associated space point, the identifiers
///        needed to match the hit against the truth or against another run of
///        the pattern recognition.
class RootMuonGlobalPatternWriter : public WriterT<MuonGlobalPatternContainer> {
 public:
  /// Configuration of the writer
  struct Config {
    /// Input global pattern collection to write.
    std::string inputPatterns{};
    /// Input space point collection the patterns were built from. It is needed
    /// to record where each pattern hit sits in the input container.
    std::string inputSpacePoints{};
    /// Path to the output file.
    std::string filePath{};
    /// Output file access mode.
    std::string fileMode{"RECREATE"};
    /// Name of the tree within the output file.
    std::string treeName{"muonGlobalPatterns"};
  };

  /// Construct the pattern writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootMuonGlobalPatternWriter(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootMuonGlobalPatternWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] patterns are the patterns to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MuonGlobalPatternContainer& patterns) override;

  Config m_cfg{};
  std::unique_ptr<TFile> m_file{};
  TTree* m_tree{};

  ReadDataHandle<MuonSpacePointContainer> m_inputSpacePoints{
      this, "InputSpacePoints"};

  mutable std::mutex m_mutex{};
  /// @brief Event identifier.
  std::uint32_t m_eventId{0};
  /// @brief Expanded sector the pattern was built in.
  ///        Stored as std::int16_t so ROOT can write the branch without a
  ///        custom dictionary, which a vector of std::int8_t would need.
  std::vector<std::int16_t> m_sector{};
  /// @brief Direction of the pattern
  std::vector<float> m_theta{};
  std::vector<float> m_phi{};
  /// @brief Counts of precision / trigger / phi layers
  std::vector<std::uint8_t> m_nPrecisionLayers{};
  std::vector<std::uint8_t> m_nTriggerLayers{};
  std::vector<std::uint8_t> m_nPhiLayers{};
  /// @brief Mean over eta hits of the square of their normalized residual
  std::vector<float> m_meanNormResidual2{};
  /// @brief Number of space points associated with the pattern
  std::vector<std::uint16_t> m_nHits{};
  /// @brief Index of the pattern the hit belongs to
  std::vector<std::uint16_t> m_hitPatternIdx{};
  /// @brief Station the hit was recorded in
  std::vector<std::uint8_t> m_hitStation{};
  /// @brief Geometry identifier of the associated surface
  std::vector<Acts::GeometryIdentifier::Value> m_hitGeometryId{};
  /// @brief Muon identifier of the hit
  std::vector<std::uint32_t> m_hitMuonId{};
  /// @brief Position of the hit in the input space point container
  std::vector<std::uint16_t> m_hitBucketId{};
  std::vector<std::uint16_t> m_hitIndexInBucket{};
};

}  // namespace ActsExamples
