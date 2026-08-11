// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/GaussianHistogramFit.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace ActsExamples {

namespace {

using Acts::Experimental::Histogram1;

/// The bin centres and contents entering a single fit, flattened out of the
/// histogram. Only bins with non-zero content are kept: ROOT's `"SQ0"` gives
/// zero-content bins zero error (Poisson counting, `GetBinError = sqrt(n)`)
/// and drops them from the least-squares sum, so keeping them here would add
/// spurious constraints that ROOT's fit does not have.
struct Bins {
  std::vector<double> centres;
  std::vector<double> counts;
};

/// Collect the non-empty bins whose centre lies within `[xMin, xMax]`
///
/// The closed interval and the use of bin centres match how ROOT restricts a
/// fit range, so that a range boundary lands on the same side of a bin here as
/// it does there.
Bins selectBins(const Histogram1& hist, double xMin, double xMax) {
  const auto& axis = hist.histogram().axis(0);

  Bins bins;
  for (int i = 0; i < axis.size(); ++i) {
    const double centre = 0.5 * (axis.bin(i).lower() + axis.bin(i).upper());
    if (centre < xMin || centre > xMax) {
      continue;
    }
    const double count = hist.binContent({i});
    if (count <= 0) {
      continue;
    }
    bins.centres.push_back(centre);
    bins.counts.push_back(count);
  }
  return bins;
}

/// Starting values for the search, following ROOT's `H1InitGaus`: the
/// count-weighted mean and RMS of the selected bins.
///
/// @return `(mean, sigma)`, or `std::nullopt` if no sensible seed exists
std::optional<Acts::Vector2> initialGuess(const Bins& bins) {
  const std::size_t n = bins.centres.size();
  const Eigen::Map<const Eigen::VectorXd> x(bins.centres.data(), n);
  const Eigen::Map<const Eigen::VectorXd> w(bins.counts.data(), n);

  const double total = w.sum();
  if (!(total > 0)) {
    return std::nullopt;
  }

  const double mean = w.dot(x) / total;
  const double variance =
      (w.array() * x.array().square()).sum() / total - mean * mean;
  // A vanishing or negative variance means the counts sit in a single bin, or
  // rounding has eaten the spread; fall back on a fraction of the fit range
  const double sigma = (variance > 0) ? std::sqrt(variance)
                                      : 0.25 * (x.maxCoeff() - x.minCoeff());
  if (!(sigma > 0) || !std::isfinite(mean)) {
    return std::nullopt;
  }

  return Acts::Vector2{mean, sigma};
}

/// Variable projection: at fixed `(mean, sigma)`, the amplitude minimising
/// the chi-square has the closed form `A = S1 / S2`,
/// `S1 = sum_i g_i`, `S2 = sum_i g_i^2 / n_i`, `g_i = exp(-z_i^2 / 2)`,
/// `z_i = (x_i - mean) / sigma`. Searching only `(mean, sigma)` with the
/// amplitude eliminated this way removes a near-degenerate direction that a
/// joint 3-parameter search can wander into (amplitude and a far-away,
/// wide Gaussian trading off along an almost-flat valley of the chi-square)
/// -- the reason this projection, not a joint fit, is used here.
///
/// @return `std::nullopt` if every included bin's model value underflows to
///         zero at this `(mean, sigma)`, i.e. `S2 == 0`
std::optional<double> profiledAmplitude(const Bins& bins, double mean,
                                        double sigma) {
  double s1 = 0;
  double s2 = 0;
  for (std::size_t i = 0; i < bins.centres.size(); ++i) {
    const double z = (bins.centres[i] - mean) / sigma;
    const double g = std::exp(-0.5 * z * z);
    s1 += g;
    s2 += g * g / bins.counts[i];
  }
  if (!(s2 > 0)) {
    return std::nullopt;
  }
  return s1 / s2;
}

/// Accumulate the Gauss-Newton normal equations `J^T J`, `J^T r` and the
/// chi-square of the full 3-parameter model at `p = (A, m, s)`, optionally
/// also accumulating the correction `S` that turns `J^T J` into half the
/// true chi-square Hessian, `J^T J - S`.
///
/// The residual is `r_i = (n_i - A g_i) / sqrt(n_i)`,
/// `g_i = exp(-z_i^2 / 2)`, `z_i = (x_i - m) / s`, with derivatives
/// `dg/dA = g_i`, `dg/dm = A g_i z_i / s`, `dg/ds = A g_i z_i^2 / s`.
///
/// Gauss-Newton (`J^T J`) drops the term coming from the model's own
/// curvature: `d^2(chi^2)/2 = J^T J - sum_i e_i r_i H_i`, `e_i = 1/sqrt(n_i)`,
/// `H_i` the Hessian of the model at bin `i`. That term vanishes only if
/// every residual is negligible; ROOT's MINUIT includes it (it differentiates
/// the actual chi-square numerically, not a linearised model of it), so
/// matching ROOT's parameter errors on anything but a near-perfect fit
/// requires it too -- but it is only needed once, for the final covariance,
/// not on every search step. `wantHessianCorrection` is false on every
/// Levenberg-Marquardt iteration below (only chi-square/JtJ/Jtr drive the
/// search) and true on the one call made afterwards for the final
/// covariance.
///
/// @return The chi-square at `p`
double normalEquations(const Bins& bins, const Acts::Vector3& p,
                       Acts::SquareMatrix3& jtj, Acts::Vector3& jtr,
                       bool wantHessianCorrection,
                       Acts::SquareMatrix3& hessianCorrection) {
  const double amplitude = p(0);
  const double mean = p(1);
  const double sigma = p(2);
  const double sigmaSq = sigma * sigma;

  jtj.setZero();
  jtr.setZero();
  hessianCorrection.setZero();
  double chiSquare = 0;

  for (std::size_t i = 0; i < bins.centres.size(); ++i) {
    const double count = bins.counts[i];
    const double invSigmaCount = 1.0 / std::sqrt(count);

    const double z = (bins.centres[i] - mean) / sigma;
    const double zSq = z * z;
    const double g = std::exp(-0.5 * zSq);
    const double model = amplitude * g;

    const Acts::Vector3 jacobianRow{
        g * invSigmaCount, amplitude * g * z / sigma * invSigmaCount,
        amplitude * g * zSq / sigma * invSigmaCount};
    const double residual = (count - model) * invSigmaCount;

    jtj += jacobianRow * jacobianRow.transpose();
    jtr += jacobianRow * residual;
    chiSquare += residual * residual;

    if (wantHessianCorrection) {
      // Second derivatives of the model A*g(m, s) w.r.t. (A, m, s)
      const double dAdm = g * z / sigma;
      const double dAds = g * zSq / sigma;
      const double dmdm = amplitude * g * (zSq - 1) / sigmaSq;
      const double dmds = amplitude * g * (z * zSq - 2 * z) / sigmaSq;
      const double dsds = amplitude * g * (zSq * zSq - 3 * zSq) / sigmaSq;

      Acts::SquareMatrix3 hessian;
      // clang-format off
      hessian << 0,    dAdm, dAds,
                dAdm, dmdm, dmds,
                dAds, dmds, dsds;
      // clang-format on
      hessianCorrection += (invSigmaCount * residual) * hessian;
    }
  }

  return chiSquare;
}

}  // namespace

std::optional<HistogramFitResult> gaussianHistogramFit(
    const Histogram1& hist, std::optional<HistogramFitRange> range) {
  constexpr std::size_t minNonEmptyBins = 3;
  constexpr std::size_t maxIterations = 50;
  constexpr double relativeTolerance = 1e-8;

  const double infinity = std::numeric_limits<double>::infinity();
  const auto [xMin, xMax] =
      range.value_or(HistogramFitRange{-infinity, infinity});
  const Bins bins = selectBins(hist, xMin, xMax);
  if (bins.centres.size() < minNonEmptyBins) {
    return std::nullopt;
  }

  const std::optional<Acts::Vector2> seed = initialGuess(bins);
  if (!seed.has_value()) {
    return std::nullopt;
  }

  // Levenberg-Marquardt on the amplitude-profiled (mean, sigma) least
  // squares. At each trial point the amplitude is re-profiled
  // (profiledAmplitude) and the full 3-parameter normal equations evaluated
  // there; the (mean, sigma) block of J^T J / J^T r is exactly the profiled
  // 2-parameter normal equations (same Jacobian columns, same residuals), so
  // no separate profiled-objective function is needed. The damped normal
  // equations `(J^T J + lambda * diag(J^T J)) delta = J^T r` interpolate
  // between a Gauss-Newton step (lambda -> 0, fast near the minimum) and a
  // small gradient-descent step (lambda large, safe far from it).
  Acts::Vector2 p = *seed;

  const auto profiledStep =
      [&bins](const Acts::Vector2& mean_sigma, Acts::SquareMatrix2& jtj,
              Acts::Vector2& jtr) -> std::optional<double> {
    const std::optional<double> amplitude =
        profiledAmplitude(bins, mean_sigma(0), mean_sigma(1));
    if (!amplitude.has_value()) {
      return std::nullopt;
    }
    Acts::SquareMatrix3 fullJtj;
    Acts::Vector3 fullJtr;
    Acts::SquareMatrix3 unusedCorrection;
    const double chiSquare = normalEquations(
        bins, {*amplitude, mean_sigma(0), mean_sigma(1)}, fullJtj, fullJtr,
        /*wantHessianCorrection=*/false, unusedCorrection);
    jtj = fullJtj.block<2, 2>(1, 1);
    jtr = fullJtr.segment<2>(1);
    return chiSquare;
  };

  Acts::SquareMatrix2 jtj;
  Acts::Vector2 jtr;
  std::optional<double> chiSquare = profiledStep(p, jtj, jtr);
  if (!chiSquare.has_value()) {
    return std::nullopt;
  }
  double lambda = 1e-3;

  bool converged = false;
  for (std::size_t iter = 0; iter < maxIterations && !converged; ++iter) {
    const Acts::SquareMatrix2 damped =
        jtj + lambda * jtj.diagonal().asDiagonal().toDenseMatrix();
    const Acts::Vector2 delta = damped.ldlt().solve(jtr);
    if (!delta.allFinite()) {
      return std::nullopt;
    }

    const Acts::Vector2 trial = p + delta;
    if (!(trial(1) > 0)) {
      // A non-positive trial sigma is never an acceptable step; treat it like
      // a failed step and increase the damping.
      lambda *= 10;
      continue;
    }

    Acts::SquareMatrix2 trialJtj;
    Acts::Vector2 trialJtr;
    const std::optional<double> trialChiSquare =
        profiledStep(trial, trialJtj, trialJtr);

    if (trialChiSquare.has_value() && *trialChiSquare < *chiSquare) {
      const double improvement = *chiSquare - *trialChiSquare;
      p = trial;
      jtj = trialJtj;
      jtr = trialJtr;
      lambda = std::max(lambda * 0.1, 1e-12);
      converged = improvement < relativeTolerance * std::max(1.0, *chiSquare);
      chiSquare = trialChiSquare;
    } else {
      lambda *= 10;
    }
  }

  const double mean = p(0);
  const double sigma = p(1);
  const std::optional<double> amplitude = profiledAmplitude(bins, mean, sigma);
  if (!p.allFinite() || !(sigma > 0) || !amplitude.has_value()) {
    return std::nullopt;
  }

  // The (mean, sigma) block of the full 3-parameter covariance is the Schur
  // complement of the amplitude row/column in half the true chi-square
  // Hessian, `J^T J - S` -- a standard identity for profiled least squares,
  // and exactly the block ROOT's `ParError` reports. Its inverse is already
  // the chi-square / MINUIT "UP = 1" covariance by construction, no extra
  // scaling needed.
  Acts::SquareMatrix3 fullJtj;
  Acts::Vector3 fullJtr;
  Acts::SquareMatrix3 hessianCorrection;
  normalEquations(bins, {*amplitude, mean, sigma}, fullJtj, fullJtr,
                  /*wantHessianCorrection=*/true, hessianCorrection);
  const Acts::SquareMatrix3 halfHessian = fullJtj - hessianCorrection;

  const double has = halfHessian(0, 0);
  if (!(has > 0)) {
    return std::nullopt;
  }
  const Acts::Vector2 ham = halfHessian.block<2, 1>(1, 0);
  const Acts::SquareMatrix2 schurComplement =
      halfHessian.block<2, 2>(1, 1) - ham * ham.transpose() / has;

  const Acts::SquareMatrix2 covariance = schurComplement.inverse();
  const double meanVariance = covariance(0, 0);
  const double sigmaVariance = covariance(1, 1);
  if (!(meanVariance > 0) || !(sigmaVariance > 0)) {
    return std::nullopt;
  }

  return HistogramFitResult{mean, sigma, std::sqrt(meanVariance),
                            std::sqrt(sigmaVariance)};
}

}  // namespace ActsExamples
