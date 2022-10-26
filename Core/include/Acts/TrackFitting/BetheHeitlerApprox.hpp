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

#include <array>
#include <fstream>
#include <mutex>
#include <random>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace Acts {

namespace detail {

struct GaussianComponent {
  ActsScalar weight, mean, var;
};

/// Transform a gaussian component to a space where all values are defined from
/// [-inf, inf]
void transformComponent(const GaussianComponent &cmp,
                        double &transformed_weight, double &transformed_mean,
                        double &transformed_var) {
  const auto &[weight, mean, var] = cmp;

  transformed_weight = std::log(weight) - std::log(1 - weight);
  transformed_mean = std::log(mean) - std::log(1 - mean);
  transformed_var = std::log(var);
}

/// Transform a gaussian component back from the [-inf, inf]-space to the usual
/// space
auto inverseTransformComponent(double transformed_weight,
                               double transformed_mean,
                               double transformed_var) {
  GaussianComponent cmp;
  cmp.weight = 1. / (1 + std::exp(-transformed_weight));
  cmp.mean = 1. / (1 + std::exp(-transformed_mean));
  cmp.var = std::exp(transformed_var);

  return cmp;
}

template <std::size_t NComponent>
class GaussianMixtureModelPDF;

template <std::size_t NComponents>
class GaussianMixtureModelCDF {
  std::array<GaussianComponent, NComponents> m_components;

 public:
  GaussianMixtureModelCDF(
      const std::array<GaussianComponent, NComponents> &components)
      : m_components(components) {}

  double operator()(double x) {
    auto cdf = [&](double mu, double sigma_squared) -> double {
      return 0.5 * (1.0 + std::erf((x - mu) / (std::sqrt(2 * sigma_squared))));
    };

    double sum = 0;

    for (const auto [weight, mean, sigma_squared] : m_components) {
      sum += weight * cdf(mean, sigma_squared);
    }

    return sum;
  }

  auto pdf() const { GaussianMixtureModelPDF<NComponents>{m_components}; }
};

/// Compute the value of the gaussian mixture at x
template <std::size_t NComponents>
class GaussianMixtureModelPDF {
  std::array<GaussianComponent, NComponents> m_components;

 public:
  GaussianMixtureModelPDF(
      const std::array<GaussianComponent, NComponents> &components)
      : m_components(components) {}

  double operator()(double x) {
    auto gaussian = [&](double mu, double sigma_squared) -> double {
      return (1. / std::sqrt(sigma_squared * 2 * M_PI)) *
             std::exp((-(x - mu) * (x - mu)) / (2 * sigma_squared));
    };

    double sum = 0;

    for (const auto [weight, mean, sigma_squared] : m_components) {
      sum += weight * gaussian(mean, sigma_squared);
    }

    return sum;
  }

  auto cdf() const { return GaussianMixtureModelCDF{m_components}; }
};

class BetheHeitlerCDF {
  double m_thickness;

 public:
  BetheHeitlerCDF(double thicknessInX0) : m_thickness(thicknessInX0) {}

  double operator()(double x) {
    if (x <= 0.0) {
      return 0.0;
    }
    if (x >= 1.0) {
      return 1.0;
    }

    auto c = m_thickness / std::log(2);
    return 1. - boost::math::gamma_p(c, -std::log(x));
  }
};

/// Compute the Bethe-Heitler Distribution on a value x
class BetheHeitlerPDF {
  double m_thickness;

 public:
  BetheHeitlerPDF(double thicknessInX0) : m_thickness(thicknessInX0) {}

  double operator()(double x) {
    if (x <= 0.0 || x >= 1.0) {
      return 0.0;
    }

    auto c = m_thickness / std::log(2);
    return std::pow(-std::log(x), c - 1) / std::tgamma(c);
  }

  BetheHeitlerCDF cdf() const { return BetheHeitlerCDF{m_thickness}; }
};

/// Integrand for the CDF distance
template <std::size_t NComponents>
class CDFIntegrant {
  GaussianMixtureModelCDF<NComponents> m_mixture;
  BetheHeitlerCDF m_distribution;

 public:
  CDFIntegrant(const GaussianMixtureModelCDF<NComponents> &mixture,
               const BetheHeitlerCDF &dist)
      : m_mixture(mixture), m_distribution(dist) {}

  double operator()(double x) {
    return std::abs(m_mixture(x) - m_distribution(x));
  }
};

template <int NComponents>
class KLIntegrant {
  GaussianMixtureModelPDF<NComponents> m_mixture;
  BetheHeitlerPDF m_distribution;

 public:
  KLIntegrant(const GaussianMixtureModelPDF<NComponents> &mixture,
              const BetheHeitlerPDF &dist)
      : m_mixture(mixture), m_distribution(dist) {}

  double operator()(double x) {
    auto dist_x = m_distribution(x) == 0. ? std::numeric_limits<double>::min()
                                          : m_distribution(x);
    return std::log(dist_x / m_mixture(x)) * dist_x;
  }
};

template <std::size_t NComponents>
class FlatCache {
  using Value =
      std::pair<double, std::array<detail::GaussianComponent, NComponents>>;
  std::vector<Value> m_cache;
  const double m_tolerance;

 public:
  FlatCache(double tol) : m_tolerance(tol) {}

  void insert(const Value &val) {
    m_cache.push_back(val);
    std::sort(m_cache.begin(), m_cache.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
  }

  std::optional<typename Value::second_type> findApprox(double x) const {
#if 1
    auto begin = std::lower_bound(
        m_cache.begin(), m_cache.end(), x - m_tolerance,
        [](const auto &p, const auto &v) { return p.first < v; });
    auto end = std::upper_bound(
        m_cache.begin(), m_cache.end(), x + m_tolerance,
        [](const auto &v, const auto &p) { return v < p.first; });

    if (begin == m_cache.end() or begin == end) {
      return std::nullopt;
    } else {
      auto ret =
          std::min_element(begin, end, [&](const auto &a, const auto &b) {
            return std::abs(a.first - m_tolerance) <
                   std::abs(b.first - m_tolerance);
          });
      return ret->second;
    }
#else
    auto cached = std::min_element(
        m_cache.begin(), m_cache.end(), [&](const auto &pa, const auto &pb) {
          return std::abs(x - pa.first) < std::abs(x - pb.first);
        });

    if (cached != m_cache.end() and std::abs(cached->first - x) < m_tolerance) {
      return cached->second;
    } else {
      return std::nullopt;
    }
#endif
  }
};
}  // namespace detail

/// This class approximates the Bethe-Heitler with only one component. The
struct BetheHeitlerApproxSingleCmp {
  /// Returns the number of components the returned mixture will have
  constexpr auto numComponents() const { return 1; }

  /// Checks if an input is valid for the parameterization. Since this is for
  /// debugging, it always returns false
  ///
  /// @param x input in terms of x/x0
  constexpr bool validXOverX0(ActsScalar) const { return false; }

  /// Returns array with length 1 containing a 1-component-representation of the
  /// Bethe-Heitler-Distribution
  static auto mixture(const ActsScalar x) {
    std::array<detail::GaussianComponent, 1> ret{};

    ret[0].weight = 1.0;

    const double c = x / std::log(2);
    ret[0].mean = std::pow(2, -c);
    ret[0].var = std::pow(3, -c) - std::pow(4, -c);

    // ret[0].mean = std::exp(-1. * x);
    // ret[0].var =
    //     std::exp(-1. * x * std::log(3.) / std::log(2.)) - std::exp(-2. * x);

    return ret;
  }
};

/// This class approximates the Bethe-Heitler distribution as a gaussian
/// mixture. To enable an approximation for continous input variables, the
/// weights, means and variances are internally parametrized as a Nth order
/// polynomial.
template <int NComponents, int PolyDegree>
class AtlasBetheHeitlerApprox {
  static_assert(NComponents > 0);
  static_assert(PolyDegree > 0);

 public:
  struct PolyData {
    std::array<ActsScalar, PolyDegree + 1> weightCoeffs, meanCoeffs, varCoeffs;
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
  /// Construct the Bethe-Heitler approximation description
  ///
  /// @param low_data data for the lower x/x0 range
  /// @param high_data data for the higher x/x0 range
  /// @param transform wether the data need to be transformed (see Atlas code)
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
  /// @param x input in terms of x/x0
  constexpr bool validXOverX0(ActsScalar x) const { return x < higherLimit; }

  /// Generates the mixture from the polynomials and reweights them, so
  /// that the sum of all weights is 1
  ///
  /// @param x The input in terms of x/x0 (pathlength in terms of radiation length)
  auto mixture(ActsScalar x) const {
    // Build a polynom
    auto poly = [](ActsScalar xx,
                   const std::array<ActsScalar, PolyDegree + 1> &coeffs) {
      ActsScalar sum{0.};
      for (const auto c : coeffs) {
        sum = xx * sum + c;
      }
      throw_assert(std::isfinite(sum), "polynom result not finite");
      return sum;
    };

    // Lambda which builds the components
    auto make_mixture = [&](const Data &data, double xx, bool transform) {
      // Value initialization should garanuee that all is initialized to zero
      std::array<detail::GaussianComponent, NComponents> ret{};
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
      std::array<detail::GaussianComponent, NComponents> ret{};

      ret[0].weight = 1.0;
      ret[0].mean = 1.0;  // p_initial = p_final
      ret[0].var = 0.0;

      return ret;
    }
    // Return single gaussian approximation
    if (x < singleGaussianLimit) {
      std::array<detail::GaussianComponent, NComponents> ret{};
      ret[0] = BetheHeitlerApproxSingleCmp::mixture(x)[0];
      return ret;
    }
    // Return a component representation for lower x0
    if (x < lowerLimit) {
      return make_mixture(m_low_data, x, m_low_transform);
    }
    // Return a component representation for higher x0
    else {
      // Cap the x because beyond the parameterization goes wild
      const auto high_x = std::min(higherLimit, x);
      return make_mixture(m_high_data, high_x, m_high_transform);
    }
  }

  /// Loads a parameterization from a file according to the Atlas file
  /// description
  ///
  /// @param low_parameters_path Path to the foo.par file that stores
  /// the parameterization for low x/x0
  /// @param high_parameters_path Path to the foo.par file that stores
  /// the parameterization for high x/x0
  static auto loadFromFile(const std::string &low_parameters_path,
                           const std::string &high_parameters_path) {
    auto read_file = [](const std::string &filepath) {
      std::ifstream file(filepath);

      if (!file) {
        throw std::invalid_argument("Could not open '" + filepath + "'");
      }

      std::size_t n_cmps, degree;
      bool transform_code;

      file >> n_cmps >> degree >> transform_code;

      if (NComponents != n_cmps) {
        throw std::invalid_argument("Wrong number of components in '" +
                                    filepath + "'");
      }

      if (PolyDegree != degree) {
        throw std::invalid_argument("Wrong wrong polynom order in '" +
                                    filepath + "'");
      }

      if (!transform_code) {
        throw std::invalid_argument("Transform-code is required in '" +
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

template <std::size_t NComponents>
struct DefaultNext {
  auto operator()(std::array<double, 3 * NComponents> ps,
                  std::mt19937 &gen) const {
    auto val_dist = std::uniform_real_distribution{-0.5, 0.5};
    for (auto &p : ps) {
      p += val_dist(gen);
    }

    return ps;
  }
};

/// This class does the approximation by minimizing the CDF distance
/// individually for each point. Probably very slow, but good vor validating.
template <std::size_t NComponents, typename next_t = DefaultNext<NComponents>>
class BetheHeitlerSimulatedAnnealingMinimizer {
  std::vector<double> m_temperatures;
  std::array<detail::GaussianComponent, NComponents> m_startValue;
  std::shared_ptr<std::mt19937> m_gen;
  next_t m_next;

  // save time by caching
  mutable detail::FlatCache<NComponents> m_cache;

 public:
  BetheHeitlerSimulatedAnnealingMinimizer(
      const std::vector<double> &temperatures,
      const std::array<detail::GaussianComponent, NComponents> &startValue,
      std::shared_ptr<std::mt19937> gen, const next_t &next = next_t{})
      : m_temperatures(temperatures),
        m_startValue(startValue),
        m_gen(gen),
        m_next(next),
        m_cache(0.001) {}

  // Make a start value that is reasonable by distributing the means with to
  // geometric series
  static auto makeStartValue() {
    std::array<detail::GaussianComponent, NComponents> mixture{};
    double m = 0.;
    for (auto i = 0ul; i < mixture.size(); ++i) {
      m += std::pow(0.5, i + 1);

      mixture[i].weight = 1. / (mixture.size() - i);
      mixture[i].mean = m;
      mixture[i].var = 0.005 / (i + 1);
    }

    detail::normalizeWeights(mixture,
                             [](auto &a) -> double & { return a.weight; });
    return mixture;
  }

  /// Returns the number of components the returned mixture will have
  constexpr auto numComponents() const { return NComponents; }

  /// Checks if an input is valid for the parameterization. Since we do the fit
  /// on-the-fly, always return true (but of course the Bethe-Heitler model does
  /// not apply to arbitrary thick materials)
  ///
  /// @param x input in terms of x/x0
  constexpr bool validXOverX0(ActsScalar) const { return true; }

  /// Performes a simulated-annealing minimization of the CDF distance at the
  /// given x/x0.
  ///
  /// @param x The input in terms of x/x0 (pathlength in terms of radiation length)
  auto mixture(const ActsScalar x,
               std::vector<double> *history = nullptr) const {
    if( history ) {
    history->reserve(m_temperatures.size());
    }

    const auto singleCmpApprox = BetheHeitlerApproxSingleCmp::mixture(x)[0];
    const double &E1 = singleCmpApprox.mean;
    const double &E2 = singleCmpApprox.var;

    if (auto cached = m_cache.findApprox(x); cached) {
      return *cached;
    }

    // Helper function to do the integration
    auto integrate = [&](const auto &mixture) {
      return boost::math::quadrature::trapezoidal(
#if 1
          detail::CDFIntegrant<NComponents>{
              detail::GaussianMixtureModelCDF<NComponents>{mixture},
              detail::BetheHeitlerCDF{x}},
#else
          detail::KLIntegrant<NComponents>{
              detail::GaussianMixtureModelPDF<NComponents>{mixture},
              detail::BetheHeitlerPDF{x}},

#endif
          -0.5, 1.5);
    };

    // Transform to coordinates defined in [-inf, inf]
    auto transform =
        [](const std::array<detail::GaussianComponent, NComponents> &cmps) {
          auto ret = std::array<double, 3 * NComponents>{};

          auto it = ret.begin();
          for (const auto &cmp : cmps) {
            detail::transformComponent(cmp, *it, *(it + 1), *(it + 2));
            it += 3;
          }

          return ret;
        };

    // Transform from coordinates defined in [-inf, inf]
    auto inv_transform = [](const std::array<double, 3 * NComponents> &cs) {
      auto ret = std::array<detail::GaussianComponent, NComponents>{};

      auto it = cs.begin();
      for (auto &cmp : ret) {
        cmp = detail::inverseTransformComponent(*it, *(it + 1), *(it + 2));
        it += 3;
      }

      return ret;
    };

    // The annealing function
    auto minimize = [&](const auto &start_value) {
      // Initialize state
      auto current_distance = integrate(start_value);
      auto current_params = transform(start_value);

      // seperately keep track of best solution
      auto best_distance = current_distance;
      auto best_params = start_value;

      for (auto T : m_temperatures) {
        if (history) {
          history->push_back(best_distance);
        }

        const auto new_params = m_next(current_params, *m_gen);
        auto trafo_params = inv_transform(new_params);
        detail::normalizeWeights(trafo_params,
                                 [](auto &a) -> double & { return a.weight; });

        std::sort(best_params.begin(), best_params.end(),
                  [](const auto &a, const auto &b) { return a.mean < b.mean; });

        const auto sum1 = std::accumulate(
            std::next(trafo_params.begin()), trafo_params.end(), 0.0,
            [](const auto &s, const auto &c) { return s + c.weight * c.mean; });

        const auto sum2 =
            std::accumulate(std::next(trafo_params.begin()), trafo_params.end(),
                            0.0, [](const auto &s, const auto &c) {
                              return s + c.weight * (c.mean * c.mean + c.var);
                            });

        trafo_params[0].mean = (1. / trafo_params[0].weight) * (E1 - sum1);

        const auto weightMeanSquared = trafo_params[0].weight *
                                       trafo_params[0].mean *
                                       trafo_params[0].mean;
        trafo_params[0].mean =
            1. / trafo_params[0].weight * (E2 - weightMeanSquared - sum2);

        const double new_distance = integrate(trafo_params);

        if (not std::isfinite(new_distance)) {
          continue;
        }

        const double p = std::exp(-(new_distance - current_distance) / T);

        if (new_distance < best_distance) {
          best_distance = new_distance;
          best_params = trafo_params;
        }

        if (new_distance < current_distance or
            p < std::uniform_real_distribution{0., 1.0}(*m_gen)) {
          current_distance = new_distance;
          current_params = new_params;
        }
      }
      return std::make_tuple(best_distance, best_params);
    };

    // Correct mean & var
#if 0
    const double E1 = std::exp(-1. * x);
    const double E2 =
        std::exp(-1. * x * std::log(3.) / std::log(2.)) - std::exp(-2. * x);
#endif

    const double threshold = 0.0025;

    auto [best_distance, best_params] = minimize(m_startValue);
#if 0
    // We want to modify the component with the smallest mean
    std::sort(best_params.begin(), best_params.end(),
              [](const auto &a, const auto &b) { return a.mean < b.mean; });

    std::cout << "best distance " << best_distance << " at x/x0 = " << x
              << std::endl;

    std::cout << "mean result: " << std::accumulate(best_params.begin(), best_params.end(), 0.0, [](const auto &s, const auto &c){ return s + c.weight * c.mean; }) << "\n\tmeans:  ";
    std::transform(best_params.begin(), best_params.end(), std::ostream_iterator<double>(std::cout, " "), [](const auto &c){ return c.mean; });
    std::cout << "\n\tweights: ";
    std::transform(best_params.begin(), best_params.end(), std::ostream_iterator<double>(std::cout, " "), [](const auto &c){ return c.weight; });
    std::cout << "\nmean true: " << singleCmpApprox.mean << "\n";

    // Evaluate these summes before modification
    const auto sum1 = std::accumulate(
        std::next(best_params.begin()), best_params.end(), 0.0,
        [](const auto &s, const auto &c) { return s + c.weight * c.mean; });

    const auto sum2 =
        std::accumulate(std::next(best_params.begin()), best_params.end(), 0.0,
                        [](const auto &s, const auto &c) {
                          return s + c.weight * (c.mean * c.mean + c.var);
                        });

    const auto new_mean = (1. / best_params[0].weight) * (E1 - sum1);

    const auto weightMeanSquared = best_params[0].weight * new_mean * new_mean;
    const auto new_var =
        1. / best_params[0].weight * (E2 - weightMeanSquared - sum2);

    auto valid = [](const detail::GaussianComponent &cmp) {
      return cmp.var > 0. && cmp.mean < 1. && cmp.var > 0.0 &&
             std::isfinite(cmp.weight);
    };

    // TODO this can not be the case if the fit is bad for very low thicknesses
    if (valid({best_params[0].weight, new_mean, new_var})) {
      best_params[0].mean = new_mean;
      best_params[0].var = new_var;
      std::cout << "successfull corrected moments\n";
    } else {
      std::cout << "could not correct moments: " << best_params[0].mean << "+-" << best_params[0].var << " -> " << new_mean << "+-" << new_var << "\n";
      auto foo = best_params;
      foo[0].mean = new_mean;
      foo[0].var = new_var;
      std::cout << "corrected mean result: " << std::accumulate(foo.begin(), foo.end(), 0.0, [](const auto &s, const auto &c){ return s + c.weight * c.mean; }) << "\n";
    }
    std::cout << "-----------------------\n";

    detail::normalizeWeights(best_params,
                             [](auto &c) -> double & { return c.weight; });

    throw_assert(std::all_of(best_params.begin(), best_params.end(), valid),
                 "bal");
    throw_assert(detail::weightsAreNormalized(
                     best_params, [](const auto &c) { return c.weight; }),
                 "not normalized");
#endif
    // Store in cache if result is good
    if (best_distance < threshold) {
      m_cache.insert({x, best_params});
    }
    return best_params;
  }
};

// template <std::size_t NComponents>
// std::map<float, std::array<GaussianComponent, NComponents>>
// BetheHeitlerSimulatedAnnealingMinimizer<NComponents>::m_cache = {};
//
// template <std::size_t NComponents>
// std::mutex
// BetheHeitlerSimulatedAnnealingMinimizer<NComponents>::m_cache_mutex = {};

/// These data are from ATLAS and allow using the GSF without loading files.
/// However, this might not be the optimal parameterization. These data come
/// this file in Athena:
/// Tracking/TrkFitter/TrkGaussianSumFilterUtils/Data/BetheHeitler_cdf_nC6_O5.par
/// These data must be transformed, so construct the AtlasBetheHeitlerApprox
/// with transforms = true
// clang-format off
constexpr static AtlasBetheHeitlerApprox<6, 5>::Data bh_cdf_cmps6_order5_data = {{
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

}  // namespace Acts
