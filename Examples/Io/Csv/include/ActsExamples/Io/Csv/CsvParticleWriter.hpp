// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <limits>
#include <string>
#include <vector>

namespace ActsExamples {

/// Write out particles in the TrackML comma-separated-value format.
///
/// This writer is restricted to outgoing particles, it is designed for
/// generated particle information.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory. Files are named
/// using the following schema
///
///     event000000001-<stem>.csv
///     event000000002-<stem>.csv
///     ...
///
/// and each line in the file corresponds to one particle.
class CsvParticleWriter final : public WriterT<SimParticleContainer> {
 public:
  struct Config {
    /// Input particles collection to write.
    std::string inputParticles;
    /// Where to place output files.
    std::string outputDir;
    /// Output filename stem.
    std::string outputStem;
    /// Number of decimal digits for floating point precision in output.
    size_t outputPrecision = std::numeric_limits<float>::max_digits10;
  };

  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  CsvParticleWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const SimParticleContainer& particles) override;

 private:
  Config m_cfg;  //!< Nested configuration struct
};

}  // namespace ActsExamples
