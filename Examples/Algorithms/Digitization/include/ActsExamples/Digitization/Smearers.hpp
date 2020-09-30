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
#include <climits>
#include <cmath>
#include <exception>
#include <random>

namespace ActsExamples {
namespace Digitization {

/// Gaussian smearing of a single parameter.
///
/// @note This smearer will smear over module boundaries
/// it has no notion of a parameter range is assumed
struct Gauss {
  std::normal_distribution<> dist{0., 1.};

  /// Construct with a @param sigma standard deviation
  Gauss(double sigma) : dist(std::normal_distribution<>(0., sigma)) {}

  /// Call operator for the SmearFunction caller interface.
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

/// Gaussian smearing of a single parameter with truncation.
///
/// In case a hit is smeared outside the range, a DigitizationError
/// indicating the truncation
struct GaussTrunc {
  std::normal_distribution<> dist{0., 1.};
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussTrunc(double sigma, const std::pair<double, double>& range_)
      : dist(std::normal_distribution<>(0., sigma)), range(range_) {}

  /// Call operator for the SmearFunction caller interface.
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
          std::pair<double, double>(svalue, dist.stddev()));
    }
    return ActsFatras::DigitizationError::SmearingOutOfRange;
  }
};

/// Gaussian smearing of a single parameter with clipping.
///
/// In case a hit is smeared outside the range, the smearing will be
/// repeated, until a maximum attempt number is reached
struct GaussClipped {
  std::normal_distribution<> dist{0., 1.};

  size_t maxAttemps = 1000;

  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussClipped(double sigma, const std::pair<double, double>& range_)
      : dist(std::normal_distribution<>(0., sigma)), range(range_) {}

  /// Call operator for the SmearFunction caller interface.
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
            std::pair<double, double>(svalue, dist.stddev()));
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

  std::uniform_real_distribution<> dist{0., 1.};

  /// Construct with a @param pitch standard deviation and @param range
  Uniform(double pitch, const std::pair<double, double>& range_)
      : binningData(Acts::open, Acts::binX,
                    (range_.second - range_.first) / pitch, range_.first,
                    range_.second) {}

  /// Constructor with a bin utility in order to get the bin borders.
  ///
  /// @param bu the bin utility which d
  Uniform(Acts::BinningData&& bd) : binningData(std::move(bd)) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call
  ///
  /// @return a Result is uniformly distributed between bin borders
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& rnd) {
    if (binningData.min < value and binningData.max > value) {
      auto bin = binningData.search(value);
      auto lower = binningData.boundaries()[bin];
      auto higher = binningData.boundaries()[bin + 1];
      double svalue = lower + (higher - lower) * dist(rnd);
      return Acts::Result<std::pair<double, double>>(
          std::pair<double, double>(svalue, (higher - lower) / std::sqrt(12.)));
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
                    (range_.second - range_.first) / pitch, range_.first,
                    range_.second) {}

  /// Constructor with a bin utility in order to get the bin borders.
  ///
  /// @param bu the bin utility within hich the parameter is allowed
  Digital(Acts::BinningData&& bd) : binningData(std::move(bd)) {}

  /// Call operator for the SmearFunction caller interface.
  ///
  /// @param value parameter to be smeared
  /// @param rnd random generator to be used for the call (unused)
  ///
  /// @return a Result is uniformly distributed between bin borders
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomEngine& /*unused*/) {
    if (binningData.min < value and binningData.max > value) {
      auto bin = binningData.search(value);
      auto lower = binningData.boundaries()[bin];
      auto higher = binningData.boundaries()[bin + 1];
      double svalue = 0.5 * (lower + higher);
      return Acts::Result<std::pair<double, double>>(
          std::pair<double, double>(svalue, (higher - lower) / std::sqrt(12.)));
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

template <size_t kSize>
using SmearFunctions =
    std::array<ActsFatras::SmearFunction<RandomEngine>, kSize>;

template <size_t kSize>
using FncTypes = std::array<int, kSize>;

template <size_t kSize>
using FncParameters = std::array<std::vector<double>, kSize>;

/// Struct to generate smearing functions from arguments.
///
struct SmearingFunctionGenerator {
  /// Generate function that unrolls the templated kSizeension.
  ///
  /// @tparam kSize The kSizeension of the smearing function array
  /// @tparam kIndex is the entry that is currently filled
  ///
  /// @param functions[in,out] The smearing functions that are generated
  /// @param types The smearing function type (Gauss as default)
  /// @param pars The parameters for the smearing function
  template <size_t kSize, size_t kIndex>
  void generateFunction(SmearFunctions<kSize>& functions,
                        const FncTypes<kSize>& types,
                        const FncParameters<kSize>& pars) noexcept(false) {
    switch (types[kIndex]) {
      case 0: {
        if (pars[kIndex].empty()) {
          throw std::invalid_argument(
              "Smearers: invalid input for Gauss smearing.");
        }
        functions[kIndex] = Gauss(pars[kIndex][0]);
      } break;

      case 1: {
        if (pars[kIndex].size() < 3) {
          throw std::invalid_argument(
              "Smearers: invalid input for truncated Gauss smearing.");
        }
        functions[kIndex] =
            GaussTrunc(pars[kIndex][0], {pars[kIndex][1], pars[kIndex][2]});
      } break;

      case 2: {
        if (pars[kIndex].size() < 3) {
          throw std::invalid_argument(
              "Smearers: invalid input for clipped Gauss smearing.");
        }
        functions[kIndex] =
            GaussClipped(pars[kIndex][0], {pars[kIndex][1], pars[kIndex][2]});
      } break;

      case 3: {
        if (pars[kIndex].size() < 3) {
          throw std::invalid_argument(
              "Smearers: invalid input for Uniform smearing.");
        }
        functions[kIndex] =
            Uniform(pars[kIndex][0], {pars[kIndex][1], pars[kIndex][2]});
      } break;

      case 4: {
        if (pars[kIndex].size() < 3) {
          throw std::invalid_argument(
              "Smearers: invalid input for Digital smearing.");
        }
        functions[kIndex] =
            Digital(pars[kIndex][0], {pars[kIndex][1], pars[kIndex][2]});
      } break;
    }

    if constexpr (kIndex > 0) {
      generateFunction<kSize, kIndex - 1>(functions, types, pars);
    }
  }

  /// Generate call for simear function generation.
  ///
  /// @tparam kSize The templated kSizeenstion type
  ///
  /// @param types The smearing function type (Gauss as default)
  /// @param pars The parameters for the smearing function
  ///
  /// @return a properly kSizeensioned resultion function array
  template <size_t kSize>
  SmearFunctions<kSize> generate(const FncTypes<kSize>& types,
                                 const FncParameters<kSize>& pars) {
    SmearFunctions<kSize> sFunctions;
    generateFunction<kSize, kSize - 1>(sFunctions, types, pars);
    return sFunctions;
  }
};

}  // namespace Digitization
}  // namespace ActsExamples
