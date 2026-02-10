// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <string>

namespace ActsExamples {

/// @class CsvSpacePointsBucketWriter
///
/// This writes one file per event containing information about the
/// buckets
///
///     event000000001-buckets.csv
///     event000000002-buckets.csv
///     ...
///
/// Intrinsically thread-safe as one file per event.
class CsvSpacePointsBucketWriter final
    : public WriterT<std::vector<SpacePointContainer>> {
 public:
  struct Config {
    /// Which bucket collection to write.
    std::string inputBuckets;
    /// Where to place output files
    std::string outputDir;
    /// Number of decimal digits for floating point precision in output.
    int outputPrecision = std::numeric_limits<float>::max_digits10;
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  CsvSpacePointsBucketWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~CsvSpacePointsBucketWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param buckets is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<SpacePointContainer>& buckets) override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
