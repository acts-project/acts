// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>

#include <boost/container/static_vector.hpp>

namespace Acts {

namespace detail {

struct GaussianComponent {
  ActsScalar weight = 0.0;
  ActsScalar mean = 0.0;
  ActsScalar var = 0.0;
};

/// Transform a gaussian component to a space where all values are defined from
/// [-inf, inf]
inline void transformComponent(const GaussianComponent &cmp,
                               double &transformed_weight,
                               double &transformed_mean,
                               double &transformed_var) {
  const auto &[weight, mean, var] = cmp;

  transformed_weight = std::log(weight) - std::log(1 - weight);
  transformed_mean = std::log(mean) - std::log(1 - mean);
  transformed_var = std::log(var);
}

/// Transform a gaussian component back from the [-inf, inf]-space to the usual
/// space
inline auto inverseTransformComponent(double transformed_weight,
                                      double transformed_mean,
                                      double transformed_var) {
  GaussianComponent cmp;
  cmp.weight = 1. / (1 + std::exp(-transformed_weight));
  cmp.mean = 1. / (1 + std::exp(-transformed_mean));
  cmp.var = std::exp(transformed_var);

  return cmp;
}

}  // namespace detail

namespace Experimental {

/// This class approximates the Bethe-Heitler with only one component. This is
/// mainly inside @ref AtlasBetheHeitlerApprox, but can also be used as the
/// only component approximation (then probably for debugging)
struct BetheHeitlerApproxSingleCmp {
  /// Returns the number of components the returned mixture will have
  constexpr auto numComponents() const { return 1; }

  /// Checks if an input is valid for the parameterization. The threshold for
  /// x/x0 is 0.002 and orientates on the values used in ATLAS
  constexpr bool validXOverX0(ActsScalar x) const {
    return x < 0.002;
    ;
  }

  /// Returns array with length 1 containing a 1-component-representation of the
  /// Bethe-Heitler-Distribution
  ///
  /// @param x pathlength in terms of the radiation length
  static auto mixture(const ActsScalar x) {
    std::array<detail::GaussianComponent, 1> ret{};

    ret[0].weight = 1.0;

    const double c = x / std::log(2);
    ret[0].mean = std::pow(2, -c);
    ret[0].var = std::pow(3, -c) - std::pow(4, -c);

    return ret;
  }
};

/// This class approximates the Bethe-Heitler distribution as a gaussian
/// mixture. To enable an approximation for continuous input variables, the
/// weights, means and variances are internally parametrized as a Nth order
/// polynomial.
template <int NComponents, int PolyDegree>
class AtlasBetheHeitlerApprox {
  static_assert(NComponents > 0);
  static_assert(PolyDegree > 0);

 public:
  struct PolyData {
    std::array<ActsScalar, PolyDegree + 1> weightCoeffs;
    std::array<ActsScalar, PolyDegree + 1> meanCoeffs;
    std::array<ActsScalar, PolyDegree + 1> varCoeffs;
  };

  using Data = std::array<PolyData, NComponents>;

  constexpr static double noChangeLimit = 0.0001;
  constexpr static double singleGaussianLimit = 0.002;
  constexpr static double lowerLimit = 0.10;
  constexpr static double higherLimit = 0.20;

 private:
  Data m_low_data;
  Data m_high_data;
  bool m_low_transform;
  bool m_high_transform;

 public:
  /// Construct the Bethe-Heitler approximation description. Additional to the
  /// coefficients of the polynomials, the information whether these values need
  /// to be transformed beforehand must be given (see ATLAS code).
  ///
  /// @param low_data data for the lower x/x0 range
  /// @param high_data data for the higher x/x0 range
  /// @param low_transform wether the low data need to be transformed
  /// @param high_transform wether the high data need to be transformed
  constexpr AtlasBetheHeitlerApprox(const Data &low_data, const Data &high_data,
                                    bool low_transform, bool high_transform)
      : m_low_data(low_data),
        m_high_data(high_data),
        m_low_transform(low_transform),
        m_high_transform(high_transform) {}

  /// Returns the number of components the returned mixture will have
  constexpr auto numComponents() const { return NComponents; }

  /// Checks if an input is valid for the parameterization
  ///
  /// @param x pathlength in terms of the radiation length
  constexpr bool validXOverX0(ActsScalar x) const { return x < higherLimit; }

  /// Generates the mixture from the polynomials and reweights them, so
  /// that the sum of all weights is 1
  ///
  /// @param x pathlength in terms of the radiation length
  auto mixture(ActsScalar x) const {
    using Array =
        boost::container::static_vector<detail::GaussianComponent, NComponents>;
    // Build a polynom
    auto poly = [](ActsScalar xx,
                   const std::array<ActsScalar, PolyDegree + 1> &coeffs) {
      ActsScalar sum{0.};
      for (const auto c : coeffs) {
        sum = xx * sum + c;
      }
      assert((std::isfinite(sum) && "polynom result not finite"));
      return sum;
    };

    // Lambda which builds the components
    auto make_mixture = [&](const Data &data, double xx, bool transform) {
      // Value initialization should garanuee that all is initialized to zero
      Array ret(NComponents);
      ActsScalar weight_sum = 0;
      for (int i = 0; i < NComponents; ++i) {
        // These transformations must be applied to the data according to ATHENA
        // (TrkGaussianSumFilter/src/GsfCombinedMaterialEffects.cxx:79)
        if (transform) {
          ret[i] = detail::inverseTransformComponent(
              poly(xx, data[i].weightCoeffs), poly(xx, data[i].meanCoeffs),
              poly(xx, data[i].varCoeffs));
        } else {
          ret[i].weight = poly(xx, data[i].weightCoeffs);
          ret[i].mean = poly(xx, data[i].meanCoeffs);
          ret[i].var = poly(xx, data[i].varCoeffs);
        }

        weight_sum += ret[i].weight;
      }

      for (int i = 0; i < NComponents; ++i) {
        ret[i].weight /= weight_sum;
      }

      return ret;
    };

    // Return no change
    if (x < noChangeLimit) {
      Array ret(1);

      ret[0].weight = 1.0;
      ret[0].mean = 1.0;  // p_initial = p_final
      ret[0].var = 0.0;

      return ret;
    }
    // Return single gaussian approximation
    if (x < singleGaussianLimit) {
      Array ret(1);
      ret[0] = BetheHeitlerApproxSingleCmp::mixture(x)[0];
      return ret;
    }
    // Return a component representation for lower x0
    if (x < lowerLimit) {
      return make_mixture(m_low_data, x, m_low_transform);
    }
    // Return a component representation for higher x0
    // Cap the x because beyond the parameterization goes wild
    const auto high_x = std::min(higherLimit, x);
    return make_mixture(m_high_data, high_x, m_high_transform);
  }

  /// Loads a parameterization from a file according to the Atlas file
  /// description
  ///
  /// @param low_parameters_path Path to the foo.par file that stores
  /// the parameterization for low x/x0
  /// @param high_parameters_path Path to the foo.par file that stores
  /// the parameterization for high x/x0
  static auto loadFromFiles(const std::string &low_parameters_path,
                            const std::string &high_parameters_path) {
    auto read_file = [](const std::string &filepath) {
      std::ifstream file(filepath);

      if (!file) {
        throw std::invalid_argument("Could not open '" + filepath + "'");
      }

      std::size_t n_cmps = 0, degree = 0;
      bool transform_code = false;

      file >> n_cmps >> degree >> transform_code;

      if (NComponents != n_cmps) {
        throw std::invalid_argument("Wrong number of components in '" +
                                    filepath + "'");
      }

      if (PolyDegree != degree) {
        throw std::invalid_argument("Wrong wrong polynom order in '" +
                                    filepath + "'");
      }

      Data data;

      for (auto &cmp : data) {
        for (auto &coeff : cmp.weightCoeffs) {
          file >> coeff;
        }
        for (auto &coeff : cmp.meanCoeffs) {
          file >> coeff;
        }
        for (auto &coeff : cmp.varCoeffs) {
          file >> coeff;
        }
      }

      return std::make_tuple(data, transform_code);
    };

    const auto [low_data, low_transform] = read_file(low_parameters_path);
    const auto [high_data, high_transform] = read_file(high_parameters_path);

    return AtlasBetheHeitlerApprox(low_data, high_data, low_transform,
                                   high_transform);
  }
};

/// Creates a @ref AtlasBetheHeitlerApprox object based on an ATLAS
/// configuration, that are stored as static data in the source code.
/// This may not be an optimal configuration, but should allow to run
/// the GSF without the need to load files
AtlasBetheHeitlerApprox<6, 5> makeDefaultBetheHeitlerApprox();

}  // namespace Experimental
}  // namespace Acts
