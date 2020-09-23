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
#include <random>

namespace ActsExamples {
namespace HitSmearers {

/// Gaussian smearing of a hit
///
/// @note This smearer will smear over module boundaries
/// it has no notion of a parameter range is assumed
struct Gauss {
  std::normal_distribution<> dist{0., 1.};

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
struct GaussTrunc {
  std::normal_distribution<> dist{0., 1.};
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussTrunc(double sigma, const std::pair<double, double>& range_)
      : dist(std::normal_distribution<>(0., sigma)), range(range_) {}

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
          std::pair<double, double>(svalue, dist.stddev()));
    }
    return ActsFatras::DigitizationError::SmearingOutOfRange;
  }
};

/// Gaussian smearing of a hit with truncation
///
/// In case a hit is smeared outside the range, the smearing will be
/// repeated, until a maximum attempt number is reached
struct GaussClipped {
  std::normal_distribution<> dist{0., 1.};
  std::pair<double, double> range = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::max()};

  /// Construct with a @param sigma standard deviation and @param range
  GaussClipped(double sigma, const std::pair<double, double>& range_)
      : dist(std::normal_distribution<>(0., sigma)), range(range_) {}

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
            std::pair<double, double>(svalue, dist.stddev()));
      }
    }
    return ActsFatras::DigitizationError::SmearingError;
  }
};

/// Uniform smearing of a hit within bounds
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

  /// Constructor with a bin utility in order to get the bin borders
  ///
  /// @param bu the bin utility which d
  Uniform(Acts::BinningData&& bd) : binningData(std::move(bd)) {}

  /// Call operator for the SmearFunction caller interface
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

template <size_t DIM>
using Functions = std::array<ActsFatras::SmearFunction<RandomEngine>, DIM>;

template <size_t DIM>
using Types = std::array<int, DIM>;

template <size_t DIM>
using Parameters = std::array<std::vector<double>, DIM>;

/// Struct to generate smearing functions from arguments
///
/// @tparam indices_t The type of free/bound indices
/// @tparam DIM The dimension of the smearing function array
struct FunctionGenerator {
  /// Template unrolling
  ///
  /// @pram ENTRY is the entry that is currently filled
  ///
  /// @param functions[in,out] The smearing functions that are generated
  /// @param types The smearing function type (Gauss as default)
  /// @param pars The parameters for the smearing function
  template <size_t DIM, size_t ENTRY>
  void generateFunction(Functions<DIM>& functions, const Types<DIM>& types,
                        const Parameters<DIM>& pars) {
    switch (types[ENTRY]) {
      case 0: {
        functions[ENTRY] = Gauss(pars[ENTRY][0]);
      } break;

      case 1: {
        functions[ENTRY] =
            GaussTrunc(pars[ENTRY][0], {pars[ENTRY][1], pars[ENTRY][2]});
      } break;

      case 2: {
        functions[ENTRY] =
            GaussClipped(pars[ENTRY][0], {pars[ENTRY][1], pars[ENTRY][2]});
      } break;

      case 3: {
        functions[ENTRY] =
            Uniform(pars[ENTRY][0], {pars[ENTRY][1], pars[ENTRY][2]});
      } break;
    }

    if constexpr (ENTRY > 0) {
      generateFunction<DIM, ENTRY - 1>(functions, types, pars);
    }
  }

  /// Generate call
  ///
  /// @param types The smearing function type (Gauss as default)
  /// @param pars The parameters for the smearing function
  ///
  /// @return a properly dimensioned resultion function array
  template <size_t DIM>
  Functions<DIM> generate(const Types<DIM>& types,
                          const Parameters<DIM>& pars) {
    Functions<DIM> sFunctions;
    generateFunction<DIM, DIM - 1>(sFunctions, types, pars);
    return sFunctions;
  }
};

}  // namespace HitSmearers
}  // namespace ActsExamples
