// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <array>
#include <fstream>

namespace Acts {
namespace detail {

inline ActsScalar logistic_sigmoid(ActsScalar x) {
  return 1. / (1 + std::exp(-x));
}

struct GaussianMixture {
  ActsScalar mean, var, weight;
};

/// This class approximates the Bethe-Heitler distribution as a gaussian
/// mixture. To enable an approximation for continous input variables, the
/// weights, means and variances are internally parametrized as a Nth order
/// polynomial.
template <int NComponents, int PolyDegree>
class BetheHeitlerApprox {
  static_assert(NComponents > 0);
  static_assert(PolyDegree > 0);

 public:
  struct PolyData {
    std::array<ActsScalar, PolyDegree + 1> weightCoeffs, meanCoeffs, varCoeffs;
  };

  using Data = std::array<PolyData, NComponents>;

 private:
  Data m_low_data;
  Data m_high_data;

  ActsScalar poly(ActsScalar x,
                  const std::array<ActsScalar, PolyDegree + 1> &coeffs) const {
    ActsScalar sum(0.);
    for (const auto c : coeffs) {
      sum = x * sum + c;
    }
    throw_assert(std::isfinite(sum), "polynom result not finite");
    return sum;
  }

 public:
  constexpr BetheHeitlerApprox(const Data &data)
      : m_low_data(data), m_high_data(data) {}
  constexpr BetheHeitlerApprox(const Data &low_data, const Data &high_data)
      : m_low_data(low_data), m_high_data(high_data) {}

  /// @brief Returns the number of components the returned mixture will have
  constexpr auto numComponents() const { return NComponents; }

  /// @brief Generates the mixture from the polynomials and reweights them, so
  /// that the sum of all weights is 1
  auto mixture(const ActsScalar x) const {
    // Value initialization should garanuee that all is initialized to zero

    // Some constants
    constexpr double singleGaussianRange = 0.0001;
    constexpr double lowerRange = 0.002;
    constexpr double higherRange = 0.10;
    constexpr double maxX0 = 0.20;

    // Lambda which builds the components
    auto make_mixture = [&](const Data &data, double xx) {
      std::array<GaussianMixture, NComponents> ret{};
      ActsScalar weight_sum = 0;
      for (int i = 0; i < NComponents; ++i) {
        // These transformations must be applied to the data according to ATHENA
        // (TrkGaussianSumFilter/src/GsfCombinedMaterialEffects.cxx:79)
        ret[i].weight = logistic_sigmoid(poly(xx, data[i].weightCoeffs));
        ret[i].mean = logistic_sigmoid(poly(xx, data[i].meanCoeffs));
        ret[i].var = std::exp(poly(xx, data[i].varCoeffs));

        weight_sum += ret[i].weight;
      }

      for (int i = 0; i < NComponents; ++i) {
        ret[i].weight /= weight_sum;
      }

      return ret;
    };

    // Return no change
    if (x < singleGaussianRange) {
      std::array<GaussianMixture, NComponents> ret{};

      ret[0].weight = 1.0;
      ret[0].mean = 1.0;  // p_initial = p_final
      ret[0].var = 0.0;

      return ret;
    }
    // Return single gaussian approximation
    if (x < lowerRange) {
      std::array<GaussianMixture, NComponents> ret{};

      ret[0].weight = 1.0;
      ret[0].mean = std::exp(-1. * x);
      ret[0].var =
          std::exp(-1. * x * std::log(3.) / std::log(2.)) - std::exp(-2. * x);

      return ret;
    }
    // Return a component representation for lower x0
    if (x < higherRange) {
      return make_mixture(m_low_data, x);
    }
    // Return a component representation for higher x0
    else {
      const auto high_x = x > maxX0 ? maxX0 : x;
      return make_mixture(m_high_data, high_x);
    }
  }
};

template <std::size_t NCmps, std::size_t Order>
auto load_bethe_heitler_data(const std::string &filepath) ->
    typename BetheHeitlerApprox<NCmps, Order>::Data {
  std::ifstream file(filepath);

  if (!file) {
    throw std::invalid_argument("Could not open '" + filepath + "'");
  }

  std::size_t n_cmps, order;
  bool transform_code;

  file >> n_cmps >> order >> transform_code;

  if (NCmps != n_cmps) {
    throw std::invalid_argument("Wrong number of components in '" + filepath +
                                "'");
  }

  if (Order != order) {
    throw std::invalid_argument("Wrong wrong polynom order in '" + filepath +
                                "'");
  }

  if (!transform_code) {
    throw std::invalid_argument("Transform-code is required in '" + filepath +
                                "'");
  }

  typename BetheHeitlerApprox<NCmps, Order>::Data data;

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

  return data;
}

/// These data are from ATLAS
/// (TrkGaussianSumFilter/Data/BetheHeitler_cdf_nC6_O5.par)
// clang-format off
constexpr static BetheHeitlerApprox<6, 5>::Data bh_cdf_cmps6_order5_data = {{
    // Component #1
    {
        {{3.74397e+004,-1.95241e+004, 3.51047e+003,-2.54377e+002, 1.81080e+001,-3.57643e+000}},
        {{3.56728e+004,-1.78603e+004, 2.81521e+003,-8.93555e+001,-1.14015e+001, 2.55769e-001}},
        {{3.73938e+004,-1.92800e+004, 3.21580e+003,-1.46203e+002,-5.65392e+000,-2.78008e+000}}
    },
    // Component #2
    {
        {{-4.14035e+004, 2.31883e+004,-4.37145e+003, 2.44289e+002, 1.13098e+001,-3.21230e+000}},
        {{-2.06936e+003, 2.65334e+003,-1.01413e+003, 1.78338e+002,-1.85556e+001, 1.91430e+000}},
        {{-5.19068e+004, 2.55327e+004,-4.22147e+003, 1.90227e+002, 9.34602e+000,-4.80961e+000}}
    },
    // Component #3
    {
        {{2.52200e+003,-4.86348e+003, 2.11942e+003,-3.84534e+002, 2.94503e+001,-2.83310e+000}},
        {{1.80405e+003,-1.93347e+003, 6.27196e+002,-4.32429e+001,-1.43533e+001, 3.58782e+000}},
        {{-4.61617e+004, 1.78221e+004,-1.95746e+003,-8.80646e+001, 3.43153e+001,-7.57830e+000}}
    },
    // Component #4
    {
        {{4.94537e+003,-2.08737e+003, 1.78089e+002, 2.29879e+001,-5.52783e+000,-1.86800e+000}},
        {{4.60220e+003,-1.62269e+003,-1.57552e+002, 2.01796e+002,-5.01636e+001, 6.47438e+000}},
        {{-9.50373e+004, 4.05517e+004,-5.62596e+003, 4.58534e+001, 6.70479e+001,-1.22430e+001}}
    },
    // Component #5
    {
        {{-1.04129e+003, 1.15222e+002,-2.70356e+001, 3.18611e+001,-7.78800e+000,-1.50242e+000}},
        {{-2.71361e+004, 2.00625e+004,-6.19444e+003, 1.10061e+003,-1.29354e+002, 1.08289e+001}},
        {{3.15252e+004,-3.31508e+004, 1.20371e+004,-2.23822e+003, 2.44396e+002,-2.09130e+001}}
    },
    // Component #6
    {
        {{1.27751e+004,-6.79813e+003, 1.24650e+003,-8.20622e+001,-2.33476e+000, 2.46459e-001}},
        {{3.64336e+005,-2.08457e+005, 4.33028e+004,-3.67825e+003, 4.22914e+001, 1.42701e+001}},
        {{-1.79298e+006, 1.01843e+006,-2.10037e+005, 1.82222e+004,-4.33573e+002,-2.72725e+001}}
    },
}};
// clang-format on

}  // namespace detail
}  // namespace Acts
