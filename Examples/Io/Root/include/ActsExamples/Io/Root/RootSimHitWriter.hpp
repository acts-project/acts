// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>

class TFile;
class TTree;

namespace ActsExamples {

/// Write out simulated hits as a flat TTree.
///
/// Each entry in the TTree corresponds to one hit for optimum writing
/// speed. The event number is part of the written data.
///
/// Safe to use from multiple writer threads. To avoid thread-saftey issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootSimHitWriter final : public WriterT<SimHitContainer> {
 public:
  struct Config {
    /// Input particle collection to write.
    std::string inputSimHits;
    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "hits";
  };

  /// Construct the particle writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootSimHitWriter(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootSimHitWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] hits are the hits to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimHitContainer& hits) final override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  /// Event identifier.
  uint32_t m_eventId;
  /// Hit surface identifier.
  uint64_t m_geometryId;
  /// Event-unique particle identifier a.k.a. barcode.
  uint64_t m_particleId;
  /// True global hit position components in mm.
  float m_tx, m_ty, m_tz;
  // True global hit time in ns.
  float m_tt;
  /// True particle four-momentum in GeV at hit position before interaction.
  float m_tpx, m_tpy, m_tpz, m_te;
  /// True change in particle four-momentum in GeV due to interactions.
  float m_deltapx, m_deltapy, m_deltapz, m_deltae;
  /// Hit index along the particle trajectory
  int32_t m_index;
  // Decoded hit surface identifier components.
  uint32_t m_volumeId;
  uint32_t m_boundaryId;
  uint32_t m_layerId;
  uint32_t m_approachId;
  uint32_t m_sensitiveId;
};

}  // namespace ActsExamples
