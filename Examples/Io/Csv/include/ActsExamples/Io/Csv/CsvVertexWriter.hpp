// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstddef>
#include <limits>
#include <string>

namespace ActsExamples {

/// Write out Vertices in a comma-separated-value format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory. Files are named
/// using the following schema
///
///     event000000001-<stem>.csv
///     event000000002-<stem>.csv
///     ...
///
/// and each line in the file corresponds to one vertex.

class CsvVertexWriter final : public WriterT<VertexContainer> {
 public:
  struct Config {
    /// Input particles collection to write.
    std::string inputVertices;
    /// Where to place output files.
    std::string outputDir;
    /// Output filename stem.
    std::string outputStem;
    /// Number of decimal digits for floating point precision in output.
    std::size_t outputPrecision = std::numeric_limits<float>::max_digits10;
  };

  /// Construct the vertex writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  CsvVertexWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] vertices are the vertices to be written
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const VertexContainer& vertices) override;

 private:
  Config m_cfg;  //!< Nested configuration struct
};

}  // namespace ActsExamples
