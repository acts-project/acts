// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <string>

namespace ActsExamples {

/// @class CsvSeedWriter
///
/// This writes one file per event containing information about the
/// seeds
///
///     event000000001-seed.csv
///     event000000002-seed.csv
///     ...
///
/// Intrinsically thread-safe as one file per event.
class CsvSeedWriter final : public WriterT<SimSeedContainer> {
 public:
  struct Config {
    /// Which seeds collection to write.
    std::string inputSeeds;
    /// Where to place output files
    std::string outputDir;
    /// Number of decimal digits for floating point precision in output.
    size_t outputPrecision = std::numeric_limits<float>::max_digits10;
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  CsvSeedWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~CsvSeedWriter() override;

  /// End-of-run hook
  ProcessCode endRun() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param seeds is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimSeedContainer& seeds) override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
