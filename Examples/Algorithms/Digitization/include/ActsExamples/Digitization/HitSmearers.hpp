// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include <climits>
#include <random>

namespace ActsExamples {
namespace HitSmearers {

/// Gaussian smearing of a hit
///
/// @note This smearer will smear over module boundaries
/// it has no notion of a parameter range is assumed
class Gauss {
  std::normal_distribution<> dist(0., 1.);

  /// Construct with a @param sigma standard deviation
  Gauss(double sigma) : dist(std::normal_distribution<>(0., sigma)) {}

  /// Call operator for the SmearFunction caller interface
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result that is always ok()
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(value + dist(rnd), dist.stddev()));
  }
};

/// Gaussian smearing of a hit with truncation
///
/// In case a hit is smeared outside the range, a DigitizationError
/// indicating the truncation
class GaussTrunc {
  std::normal_distribution<> dist(0., 1.);
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Call operator for the SmearFunction caller interface
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result that is ok() when inside range, other DigitizationError
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) {
    double svalue = value + dist(rnd);
    if (svalue >= range.first and svalue <= range.second) {
      return Acts::Result<std::pair<double, double>>(
          std::make_pair<double, double>(svalue, dist.stddev()));
    }
    return ActsFatras::DigitizationError::SmearingOutOfRange;
  }
};

/// Gaussian smearing of a hit with truncation
///
/// In case a hit is smeared outside the range, the smearing will be 
/// repeated, until a maximum attempt number is reached
class GaussClipped {
  std::normal_distribution<> dist(0., 1.);
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  size_t maxAttemps = 1000;

  /// Call operator for the SmearFunction caller interface
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @note it will smear until inside range, unless maxAttempts is reached
  ///
  /// @return a Result that is ok() when inside range, other DigitizationError
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) {
    for (size_t attempt = 0; attempt < maxAttemps; ++attempt) {
      double svalue = value + dist(rnd);
      if (svalue >= range.first and svalue <= range.second) {
        return Acts::Result<std::pair<double, double>>(
            std::make_pair<double, double>(svalue, dist.stddev()));
      }
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

/// Uniform smearing of a hit within bounds
///
/// It estimates the bin borders and smears uniformly between them
class Uniform {

    Acts::BinUtility binUtility;

    /// Constructor with a bin utility in order to get the bin borders
    ///
    /// @param bu the bin utility which d
    Uniform(Acts::BinUtility&& bu) : binUtility(std::move(bu)) {}


};

}  // namespace HitSmearers
}  // namespace ActsExamples
