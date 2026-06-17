// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <string>
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

  /// Loads a parameterization from a file according to the Atlas file
  /// description
  ///
  /// @param low_parameters_path Path to the foo.par file that stores
  /// the parameterization for low x/x0
  /// @param high_parameters_path Path to the foo.par file that stores
  /// the parameterization for high x/x0
  /// @param lowLimit the upper limit for the low x/x0-data
  /// @param highLimit the upper limit for the high x/x0-data
  /// @param clampToRange forwarded to constructor
  /// @param noChangeLimit forwarded to constructor
  /// @param singleGaussianLimit forwarded to constructor
  /// @return PolynomialBetheHeitlerApprox instance loaded from parameter files
  /// @deprecated Use loadBetheHeitlerApproxFromJson instead
  [[deprecated(
      "loadFromFiles is deprecated. Use loadBetheHeitlerApproxFromJson "
      "instead.")]]
  static PolynomialBetheHeitlerApprox loadFromFiles(
      const std::string &low_parameters_path,
      const std::string &high_parameters_path, double lowLimit,
      double highLimit, bool clampToRange, double noChangeLimit,
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
  PolynomialBetheHeitlerApprox(const Data &lowData, const Data &highData,
                               bool lowTransform, bool highTransform,
                               double lowLimit, double highLimit,
                               bool clampToRange, double noChangeLimit,
                               double singleGaussianLimit)
      : m_lowData(lowData),
        m_highData(highData),
        m_lowTransform(lowTransform),
        m_highTransform(highTransform),
        m_lowLimit(lowLimit),
        m_highLimit(highLimit),
        m_clampToRange(clampToRange),
        m_noChangeLimit(noChangeLimit),
        m_singleGaussianLimit(singleGaussianLimit) {}

  /// Returns the number of components the returned mixture will have
  /// @return Number of components in the mixture
  std::size_t maxComponents() const override {
    return std::max(m_lowData.size(), m_highData.size());
  }

  /// Checks if an input is valid for the parameterization
  ///
  /// @param xOverX0 pathlength in terms of the radiation length
  /// @return True if x/x0 is within valid range for this parameterization
  bool validXOverX0(const double xOverX0) const override {
    if (m_clampToRange) {
      return true;
    } else {
      return xOverX0 < m_highLimit;
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
  Data m_lowData;
  Data m_highData;
  bool m_lowTransform = false;
  bool m_highTransform = false;
  double m_lowLimit = 0;
  double m_highLimit = 0;
  bool m_clampToRange = false;
  double m_noChangeLimit = 0;
  double m_singleGaussianLimit = 0;
};

/// @deprecated Use PolynomialBetheHeitlerApprox instead
using AtlasBetheHeitlerApprox = PolynomialBetheHeitlerApprox;

/// @deprecated Use PolynomialBetheHeitlerApprox instead
/// @note This function is deprecated for documentation purposes. It still
///       returns PolynomialBetheHeitlerApprox which is the new default.
[[deprecated(
    "AtlasBetheHeitlerApprox is deprecated. Use PolynomialBetheHeitlerApprox "
    "instead.")]]
inline PolynomialBetheHeitlerApprox makeDefaultBetheHeitlerApprox(
    bool clampToRange = false) {
  // Tracking/TrkFitter/TrkGaussianSumFilterUtils/Data/BetheHeitler_cdf_nC6_O5.par
  // clang-format off
  static PolynomialBetheHeitlerApprox::Data cdf_cmps6_order5_data = {{
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

  return {cdf_cmps6_order5_data,
          cdf_cmps6_order5_data,
          true,
          true,
          0.2,
          0.2,
          clampToRange,
          0.0001,
          0.002};
}

/// @}

}  // namespace Acts
