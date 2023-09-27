// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <utility>

namespace ActsExamples {
namespace Digitization {

/// Exact smearing of a single parameter.
///
struct Exact {
  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result that is always ok(), and just returns
  /// the value and a stddev of 0.0
  Acts::Result<std::pair<double, double>> operator()(
      double value, RandomEngine& /*rnd*/) const {
    return std::pair{value, 0.0};
  }
};

/// Gaussian smearing of a single parameter.
///
/// @note This smearer will smear over module boundaries
/// it has no notion of a parameter range is assumed
struct Gauss {
  double sigma;

  /// Construct with a @param sigma standard deviation
  Gauss(double sigma_) : sigma{sigma_} {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result that is always ok()
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) const {
    std::normal_distribution<> dist{0, sigma};
    return std::pair{value + dist(rnd), dist.stddev()};
  }
};

/// Gaussian smearing of a single parameter with truncation.
///
/// In case a hit is smeared outside the range, a DigitizationError
/// indicating the truncation
struct GaussTrunc {
  double sigma;
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussTrunc(double sigma_, const std::pair<double, double>& range_)
      : sigma{sigma_}, range(range_) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result that is ok() when inside range, other DigitizationError
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) const {
    std::normal_distribution<> dist{0., sigma};
    double svalue = value + dist(rnd);
    if (svalue >= range.first and svalue <= range.second) {
      return std::pair{svalue, dist.stddev()};
    }
    return ActsFatras::DigitizationError::SmearingOutOfRange;
  }
};

/// Gaussian smearing of a single parameter with clipping.
///
/// In case a hit is smeared outside the range, the smearing will be
/// repeated, until a maximum attempt number is reached
struct GaussClipped {
  double sigma;

  size_t maxAttemps = 1000;

  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussClipped(double sigma_, const std::pair<double, double>& range_)
      : sigma{sigma_}, range(range_) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @note it will smear until inside range, unless maxAttempts is reached
  ///
  /// @return a Result that is ok() when inside range, other DigitizationError
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) const {
    std::normal_distribution<> dist{0., sigma};
    for (size_t attempt = 0; attempt < maxAttemps; ++attempt) {
      double svalue = value + dist(rnd);
      if (svalue >= range.first and svalue <= range.second) {
        return std::pair{svalue, dist.stddev()};
      }
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

/// Uniform smearing of a single parameter within bounds.
///
/// It estimates the bin borders and smears uniformly between them
struct Uniform {
  Acts::BinningData binningData;

  /// Construct with a @param pitch standard deviation and @param range
  Uniform(double pitch, const std::pair<double, double>& range_)
      : binningData(Acts::open, Acts::binX,
                    static_cast<size_t>((range_.second - range_.first) / pitch),
                    range_.first, range_.second) {}

  /// Constructor with a binning data in order to get the bin borders.
  ///
  /// @param bu the binning data
  Uniform(Acts::BinningData&& bd) : binningData(bd) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result is uniformly distributed between bin borders
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) const {
    if (binningData.min < value and binningData.max > value) {
      auto bin = binningData.search(value);
      auto lower = binningData.boundaries()[bin];
      auto higher = binningData.boundaries()[bin + 1];
      std::uniform_real_distribution<> dist{0., 1.};
      double svalue = lower + (higher - lower) * dist(rnd);
      return std::pair{svalue, (higher - lower) / std::sqrt(12.)};
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

/// Digital emulation of a single parameter.
///
/// It estimates the bin and gives the bin center value
struct Digital {
  Acts::BinningData binningData;

  /// Construct with a @param pitch standard deviation and @param range
  Digital(double pitch, const std::pair<double, double>& range_)
      : binningData(Acts::open, Acts::binX,
                    static_cast<size_t>((range_.second - range_.first) / pitch),
                    range_.first, range_.second) {}

  /// Constructor with a bin utility in order to get the bin borders.
  ///
  /// @param bu the bin utility within hich the parameter is allowed
  Digital(Acts::BinningData&& bd) : binningData(bd) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call (unused)
  ///
  /// @return a Result is uniformly distributed between bin borders
  Acts::Result<std::pair<double, double>> operator()(
      double value, RandomEngine& /*rnd*/) const {
    if (binningData.min < value and binningData.max > value) {
      auto bin = binningData.search(value);
      auto lower = binningData.boundaries()[bin];
      auto higher = binningData.boundaries()[bin + 1];
      double svalue = 0.5 * (lower + higher);
      return std::pair{svalue, (higher - lower) / std::sqrt(12.)};
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

}  // namespace Digitization
}  // namespace ActsExamples
