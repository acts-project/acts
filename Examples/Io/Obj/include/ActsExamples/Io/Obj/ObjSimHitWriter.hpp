// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>

namespace ActsExamples {
struct AlgorithmContext;

/// Write out a simhit collection before detector digitization
/// as wavefront obj files.
///
/// This writes one file per event containing information about the
/// global space points, momenta (before and after interaction) and hit index
/// into the configured output directory. By default it writes to the
/// current working directory. Files are named using the following schema
///
///     event000000001-<stem>.obj
///     event000000002-<stem>.obj
///     ...class ObjSimHitWriter final : public WriterT<SimHitContainer> {
class ObjSimHitWriter : public WriterT<SimHitContainer> {
 public:
  struct Config {
    /// Input sim hit collection to write.
    std::string inputSimHits;
    /// Where to place output files
    std::string outputDir;
    /// Output filename stem.
    std::string outputStem = "simhits";
    /// Number of decimal digits for floating point precision in output.
    std::size_t outputPrecision = std::numeric_limits<float>::max_digits10;
    /// Draw line connections between hits
    bool drawConnections = true;
  };

  /// Construct the particle writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  ObjSimHitWriter(const Config& config, Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~ObjSimHitWriter() override = default;

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
                     const SimHitContainer& hits) override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;
};

}  // namespace ActsExamples
