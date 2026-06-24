// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/RangeXD.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

namespace Acts {

namespace detail {

struct GaussianComponent {
  double weight = 0;
  double mean = 0;
  double var = 0;
};

/// Transform a gaussian component to a space where all values are defined from
/// [-inf, inf]
inline GaussianComponent transformComponent(const GaussianComponent &cmp) {
  GaussianComponent transformed_cmp;
  transformed_cmp.weight = std::log(cmp.weight) - std::log(1 - cmp.weight);
  transformed_cmp.mean = std::log(cmp.mean) - std::log(1 - cmp.mean);
  transformed_cmp.var = std::log(cmp.var);
  return transformed_cmp;
}

/// Transform a gaussian component back from the [-inf, inf]-space to the usual
/// space
inline GaussianComponent inverseTransformComponent(
    const GaussianComponent &transformed_cmp) {
  GaussianComponent cmp;
  cmp.weight = 1 / (1 + std::exp(-transformed_cmp.weight));
  cmp.mean = 1 / (1 + std::exp(-transformed_cmp.mean));
  cmp.var = std::exp(transformed_cmp.var);
  return cmp;
}

}  // namespace detail

/// @addtogroup track_fitting
/// @{

/// Interface for Bethe-Heitler Gaussian mixture approximations.
/// @ingroup material
class BetheHeitlerApprox {
 public:
  /// Type alias for Gaussian mixture component
  using Component = detail::GaussianComponent;

  virtual ~BetheHeitlerApprox() = default;

  /// Maximum number of components in the mixture
  /// @return Maximum number of components
  virtual std::size_t maxComponents() const = 0;

  /// Check if x/X0 value is valid for this approximation
  /// @param xOverX0 Material thickness in radiation lengths
  /// @return True if value is valid
  virtual bool validXOverX0(double xOverX0) const = 0;

  /// Compute mixture for given x/X0
  /// @param xOverX0 Material thickness in radiation lengths
  /// @param mixture Output span for mixture components
  /// @return Span of computed mixture components
  virtual std::span<Component> mixture(double xOverX0,
                                       std::span<Component> mixture) const = 0;
};

/// This class approximates the Bethe-Heitler with only one component. This is
/// mainly inside @ref PolynomialBetheHeitlerApprox, but can also be used as the
/// only component approximation (then probably for debugging)
/// @ingroup material
class BetheHeitlerApproxSingleCmp final : public BetheHeitlerApprox {
 public:
  /// Returns the number of components the returned mixture will have
  /// @return Number of components (always 1 for single component approximation)
  std::size_t maxComponents() const override { return 1; }

  /// Checks if an input is valid for the parameterization. The threshold for
  /// x/x0 is 0.002 and orientates on the values used in ATLAS
  /// @param xOverX0 The x/x0 value to check
  /// @return True if x/x0 is below the threshold for single component approximation
  bool validXOverX0(const double xOverX0) const override {
    return xOverX0 < 0.002;
  }

  /// Returns array with length 1 containing a 1-component-representation of the
  /// Bethe-Heitler-Distribution
  ///
  /// @param xOverX0 pathlength in terms of the radiation length
  /// @param mixture preallocated array to store the result
  /// @return the unmodified input span containing the single component
  std::span<Component> mixture(
      const double xOverX0, const std::span<Component> mixture) const override {
    mixture[0].weight = 1.0;

    const double c = xOverX0 / std::numbers::ln2;
    mixture[0].mean = std::pow(2, -c);
    mixture[0].var = std::pow(3, -c) - std::pow(4, -c);

    return mixture;
  }
};

/// This class approximates the Bethe-Heitler distribution as a gaussian
/// mixture. To enable an approximation for continuous input variables, the
/// weights, means and variances are internally parametrized as a Nth order
/// polynomial.
/// @ingroup material
class PolynomialBetheHeitlerApprox : public BetheHeitlerApprox {
 public:
  /// Polynomial coefficient sets for a Gaussian mixture component.
  struct PolyData {
    /// Polynomial coefficients for component weight
    std::vector<double> weightCoeffs;
    /// Polynomial coefficients for component mean
    std::vector<double> meanCoeffs;
    /// Polynomial coefficients for component variance
    std::vector<double> varCoeffs;
  };

  /// Type alias for array of polynomial data for all components
  using Data = std::vector<PolyData>;

  /// Single x/x0 range with its data and transformation flag
  struct RangeData {
    Range1D<double> range;
    Data data;
    bool transform = false;

    RangeData() = default;
    RangeData(Range1D<double> r, Data d, bool t)
        : range(r), data(std::move(d)), transform(t) {}
    RangeData(double l, double h, Data d, bool t)
        : range(l, h), data(std::move(d)), transform(t) {}
  };

  /// Construct the Bethe-Heitler approximation description with N ranges.
  /// Each range has its own data and transformation flag.
  /// The ranges will be sorted by minimum value and validated for
  /// non-overlapping.
  ///
  /// @param ranges Vector of range data
  /// @param clampToRange whether to clamp the input x/x0 to the allowed range
  /// @param noChangeLimit limit below which no change is applied
  /// @param singleGaussianLimit limit below which a single Gaussian is used
  PolynomialBetheHeitlerApprox(std::vector<RangeData> ranges, bool clampToRange,
                               double noChangeLimit,
                               double singleGaussianLimit);

  /// Construct the Bethe-Heitler approximation description with two
  /// parameterizations, one for lower ranges, one for higher ranges.
  /// Is it assumed that the lower limit of the high-x/x0 data is equal
  /// to the upper limit of the low-x/x0 data.
  ///
  /// @param lowData data for the lower x/x0 range
  /// @param highData data for the higher x/x0 range
  /// @param lowTransform whether the low data need to be transformed
  /// @param highTransform whether the high data need to be transformed
  /// @param lowLimit the upper limit for the low data
  /// @param highLimit the upper limit for the high data
  /// @param clampToRange whether to clamp the input x/x0 to the allowed range
  /// @param noChangeLimit limit below which no change is applied
  /// @param singleGaussianLimit limit below which a single Gaussian is used
  [[deprecated("Use constructor taking std::vector<RangeData> instead")]]
  PolynomialBetheHeitlerApprox(const Data &lowData, const Data &highData,
                               bool lowTransform, bool highTransform,
                               double lowLimit, double highLimit,
                               bool clampToRange, double noChangeLimit,
                               double singleGaussianLimit)
      : PolynomialBetheHeitlerApprox(
            {{0.0, lowLimit, lowData, lowTransform},
             {lowLimit, highLimit, highData, highTransform}},
            clampToRange, noChangeLimit, singleGaussianLimit) {}

  /// Returns the number of components the returned mixture will have
  /// @return Number of components in the mixture
  std::size_t maxComponents() const override {
    return std::ranges::max(
        m_ranges, {}, [](const auto &range) { return range.data.size(); });
  }

  /// Checks if an input is valid for the parameterization
  ///
  /// @param xOverX0 pathlength in terms of the radiation length
  /// @return True if x/x0 is within valid range for this parameterization
  bool validXOverX0(const double xOverX0) const override {
    if (m_clampToRange) {
      return true;
    } else {
      return xOverX0 < rangeMax();
    }
  }

  /// Generates the mixture from the polynomials and reweights them, so
  /// that the sum of all weights is 1
  ///
  /// @param xOverX0 pathlength in terms of the radiation length
  /// @param mixture preallocated array to store the result
  /// @return the potentially modified input span containing the mixture
  std::span<Component> mixture(
      double xOverX0, const std::span<Component> mixture) const override;

 private:
  double rangeMax() const { return m_ranges.back().range.max(); }

  std::vector<RangeData> m_ranges;
  bool m_clampToRange = false;
  double m_noChangeLimit = 0;
  double m_singleGaussianLimit = 0;
};

/// @deprecated Use PolynomialBetheHeitlerApprox instead
using AtlasBetheHeitlerApprox = PolynomialBetheHeitlerApprox;

/// Creates a @ref PolynomialBetheHeitlerApprox object based on a default
/// configuration, stored as static data in the source code.
/// This may not be an optimal configuration, but allows running
/// the GSF without loading external files.
/// @param clampToRange Whether to clamp values to the valid range
/// @return PolynomialBetheHeitlerApprox with default configuration parameters
/// @ingroup material
PolynomialBetheHeitlerApprox makeDefaultBetheHeitlerApprox(
    bool clampToRange = false);

/// @}

}  // namespace Acts
