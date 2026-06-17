// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <tuple>

namespace Acts {

PolynomialBetheHeitlerApprox PolynomialBetheHeitlerApprox::loadFromFiles(
    const std::string &low_parameters_path,
    const std::string &high_parameters_path, double lowLimit, double highLimit,
    bool clampToRange, double noChangeLimit, double singleGaussianLimit) {
  const auto read_file = [](const std::string &filepath) {
    std::ifstream file(filepath);

    if (!file) {
      throw std::invalid_argument("Could not open '" + filepath + "'");
    }

    std::size_t n_cmps = 0;
    std::size_t degree = 0;
    bool transform_code = false;

    file >> n_cmps >> degree >> transform_code;

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

  const auto [lowData, lowTransform] = read_file(low_parameters_path);
  const auto [highData, highTransform] = read_file(high_parameters_path);

  return {lowData,   highData,     lowTransform,  highTransform,      lowLimit,
          highLimit, clampToRange, noChangeLimit, singleGaussianLimit};
}

std::span<PolynomialBetheHeitlerApprox::Component>
PolynomialBetheHeitlerApprox::mixture(
    double xOverX0, const std::span<Component> mixture) const {
  if (m_clampToRange) {
    xOverX0 = std::clamp(xOverX0, 0.0, m_highLimit);
  }

  // Evaluate polynomial at x
  const auto poly = [](const double xx, const std::span<const double> coeffs) {
    double sum{0.};
    for (const auto c : coeffs) {
      sum = xx * sum + c;
    }
    return sum;
  };

  // Lambda which builds the components
  const auto make_mixture = [&](const Data &data, double xx,
                                bool transform) -> std::span<Component> {
    // Value initialization should garanuee that all is initialized to zero
    double weight_sum = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
      // These transformations must be applied to the data according to ATHENA
      // (TrkGaussianSumFilter/src/GsfCombinedMaterialEffects.cxx:79)
      if (transform) {
        mixture[i] = detail::inverseTransformComponent(
            {poly(xx, data[i].weightCoeffs), poly(xx, data[i].meanCoeffs),
             poly(xx, data[i].varCoeffs)});
      } else {
        mixture[i].weight = poly(xx, data[i].weightCoeffs);
        mixture[i].mean = poly(xx, data[i].meanCoeffs);
        mixture[i].var = poly(xx, data[i].varCoeffs);
      }

      weight_sum += mixture[i].weight;
    }

    for (std::size_t i = 0; i < data.size(); ++i) {
      mixture[i].weight /= weight_sum;
    }

    return {mixture.data(), data.size()};
  };

  // Return no change
  if (xOverX0 < m_noChangeLimit) {
    mixture[0].weight = 1.0;
    mixture[0].mean = 1.0;  // p_initial = p_final
    mixture[0].var = 0.0;

    return {mixture.data(), 1};
  }

  // Return single gaussian approximation
  if (xOverX0 < m_singleGaussianLimit) {
    BetheHeitlerApproxSingleCmp().mixture(xOverX0, mixture);
    return {mixture.data(), 1};
  }

  // Return a component representation for lower x0
  if (xOverX0 < m_lowLimit) {
    return make_mixture(m_lowData, xOverX0, m_lowTransform);
  }

  // Return a component representation for higher x0
  // Cap the x because beyond the parameterization goes wild
  const auto high_x = std::min(m_highLimit, xOverX0);
  return make_mixture(m_highData, high_x, m_highTransform);
}

}  // namespace Acts
