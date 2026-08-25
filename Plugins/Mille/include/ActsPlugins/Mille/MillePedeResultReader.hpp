// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <filesystem>
#include <memory>
#include <vector>

namespace ActsPlugins {

class MillePedeResultReader {
 public:
  /// @brief alignment fit result for one parameter.
  /// Packages the raw output, agnostic of the mapping
  /// to surfaces or transforms. The latter translation
  /// is expected to be done in a follow-up step.
  struct ParameterResult {
    int label = 0;     /// label of the parameter
    double val = 0;    /// fitted value
    double start = 0;  /// starting value
    double delta = 0;  /// parameter change in this iteration
    double sigma = 0;  /// uncertainty
    int nRecords = 0;  /// number of measurements affected by the parameter
  };

  /// @brief constructor - nothing to do as minimal internal state carried
  explicit MillePedeResultReader(
      std::unique_ptr<const Acts::Logger> _logger =
          Acts::getDefaultLogger("MillePedeResultReader", Acts::Logging::INFO))
      : m_logger(std::move(_logger)) {}

  /// @brief Read the alignment parameters from a result file.
  /// @param mpFile: Location of the file to parse
  Acts::Result<std::vector<ParameterResult>> readParameters(
      const std::filesystem::path& mpFile) const;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsPlugins
