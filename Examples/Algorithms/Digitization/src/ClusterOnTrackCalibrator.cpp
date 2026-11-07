// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/ClusterOnTrackCalibrator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

namespace {

/// The projection of a cluster onto one axis of the segmentation.
struct ClusterProjection {
  /// The edges of the first and the last bin of the cluster on this axis
  double firstLower = 0;
  double firstUpper = 0;
  double lastLower = 0;
  double lastUpper = 0;
  /// The summed activation of the cells in the first and the last bin
  double firstCharge = 0;
  double lastCharge = 0;
  /// The number of cells in the first and the last bin
  std::size_t firstCells = 0;
  std::size_t lastCells = 0;
  /// The widths of the bins before the first and after the last bin. Zero if
  /// the sensor ends there.
  double beforeWidth = 0;
  double afterWidth = 0;
};

/// The footprint of the predicted track in the sensor, along one axis.
struct ClusterFootprint {
  /// The extent of the path through the sensor on this axis
  double length = 0;
  /// The error of the extent, from the error of the predicted direction
  double lengthError = 0;
  /// The full length of the path through the sensor
  double pathLength = 0;
};

/// The mean and variance of the track position on one axis.
struct PosteriorMoments {
  double mean = 0;
  double variance = 0;
};

constexpr double kMinusInfinity = -std::numeric_limits<double>::infinity();

/// The log of the standard normal CDF, without underflow in the tails.
/// Both branches get clamped arguments, because the compiler can evaluate
/// both.
double logNormalCdf(double x) {
  constexpr double kTail = 8.;
  const double inner = std::clamp(x, -kTail, kTail);
  const double tail = std::min(x, -kTail);
  const double tail2 = tail * tail;
  // Asymptotic expansion of the Mills ratio
  const double asymptotic = -0.5 * tail2 - std::log(-tail) -
                            0.5 * std::log(2. * std::numbers::pi) +
                            std::log1p(-1. / tail2 + 3. / (tail2 * tail2));
  const double central =
      std::log(0.5 * std::erfc(-inner / std::numbers::sqrt2));
  return x > kTail ? 0. : (x < -kTail ? asymptotic : central);
}

/// The fraction of the segment with midpoint `mid` and length `length` in
/// the bin [lower, upper). The division is safe for a zero length, because
/// the compiler can evaluate it before the length check.
double overlapFraction(double mid, double length, double lower, double upper) {
  const double overlap =
      std::min(upper, mid + 0.5 * length) - std::max(lower, mid - 0.5 * length);
  const double fraction = std::max(overlap, 0.) /
                          std::max(length, std::numeric_limits<double>::min());
  const double inside = (lower <= mid && mid < upper) ? 1. : 0.;
  return length > 0 ? fraction : inside;
}

/// The log likelihood of the cluster for the track midpoint `mid`.
struct ClusterLikelihood {
  const ClusterProjection& projection;
  double pathLength;
  double threshold;
  double noise;
  bool analogue;

  /// The log probability that a cell with this charge fires or not
  double logFire(double charge, double sigma, bool fires) const {
    const double z = (charge - threshold) / sigma;
    return logNormalCdf(fires ? z : -z);
  }

  /// The log likelihood of the observed charge of a fired bin
  double logFired(double charge, double observed, std::size_t cells) const {
    const double sigma = noise * std::sqrt(static_cast<double>(cells));
    if (!analogue) {
      return logFire(charge, sigma, true);
    }
    const double z = (observed - charge) / sigma;
    return -0.5 * z * z;
  }

  double operator()(double mid, double length) const {
    const auto& p = projection;
    const double first =
        pathLength * overlapFraction(mid, length, p.firstLower, p.firstUpper);
    const double last =
        pathLength * overlapFraction(mid, length, p.lastLower, p.lastUpper);
    // The segment touches every fired bin
    if (first <= 0 || last <= 0) {
      return kMinusInfinity;
    }
    double result = logFired(first, p.firstCharge, p.firstCells);
    if (p.firstLower != p.lastLower) {
      result += logFired(last, p.lastCharge, p.lastCells);
    }
    // A touched neighbour stays below the threshold, an untouched one never
    // fires
    if (p.beforeWidth > 0) {
      const double before =
          pathLength * overlapFraction(mid, length,
                                       p.firstLower - p.beforeWidth,
                                       p.firstLower);
      if (before > 0) {
        result += logFire(before, noise, false);
      }
    }
    if (p.afterWidth > 0) {
      const double after =
          pathLength *
          overlapFraction(mid, length, p.lastUpper, p.lastUpper + p.afterWidth);
      if (after > 0) {
        result += logFire(after, noise, false);
      }
    }
    return result;
  }
};

/// A footprint length with its weight
using WeightedLength = std::pair<double, double>;

/// log(exp(a) + exp(b)), where either can be minus infinity
double logAddExp(double a, double b) {
  const double larger = std::max(a, b);
  const double smaller = std::min(a, b);
  return smaller == kMinusInfinity
             ? larger
             : larger + std::log1p(Acts::safeExp(smaller - larger));
}

/// Integrate the posterior of the track position on a grid.
///
/// The likelihood is the weighted sum over the footprint lengths. The grid
/// covers the positions where the segment touches the first and the last bin
/// of the cluster. If the posterior is narrow compared to the grid step, the
/// grid is refined around the peak.
///
/// @param likelihood The log likelihood of the cluster
/// @param projection The cluster on this axis
/// @param lengths The footprint lengths with their weights
///
/// @return The moments, or nothing if no position gives the cluster
std::optional<PosteriorMoments> integratePosterior(
    const ClusterLikelihood& likelihood, const ClusterProjection& projection,
    std::span<const WeightedLength> lengths) {
  constexpr std::size_t kPoints = 64;
  constexpr std::size_t kMaxPasses = 4;

  double maxLength = 0;
  for (const auto& [length, weight] : lengths) {
    maxLength = std::max(maxLength, length);
  }
  double low = projection.lastLower - 0.5 * maxLength;
  double high = projection.firstUpper + 0.5 * maxLength;
  // The hit is on the sensor, even when its footprint is clipped at an edge.
  if (projection.beforeWidth == 0) {
    low = std::max(low, projection.firstLower);
  }
  if (projection.afterWidth == 0) {
    high = std::min(high, projection.lastUpper);
  }
  if (!(low < high)) {
    return std::nullopt;
  }

  std::array<double, kPoints> logLikelihood{};
  PosteriorMoments moments;
  for (std::size_t pass = 0; pass < kMaxPasses; ++pass) {
    const double step = (high - low) / kPoints;
    const auto position = [&](std::size_t i) {
      return low + (static_cast<double>(i) + 0.5) * step;
    };

    double maxLog = kMinusInfinity;
    for (std::size_t i = 0; i < kPoints; ++i) {
      double sum = kMinusInfinity;
      for (const auto& [length, weight] : lengths) {
        sum =
            logAddExp(sum, likelihood(position(i), length) + std::log(weight));
      }
      logLikelihood[i] = sum;
      maxLog = std::max(maxLog, sum);
    }
    if (maxLog == kMinusInfinity) {
      return std::nullopt;
    }

    double norm = 0;
    double sum1 = 0;
    double sum2 = 0;
    const double origin = 0.5 * (low + high);
    for (std::size_t i = 0; i < kPoints; ++i) {
      const double w = Acts::safeExp(logLikelihood[i] - maxLog);
      const double delta = position(i) - origin;
      norm += w;
      sum1 += w * delta;
      sum2 += w * delta * delta;
    }
    const double mean = sum1 / norm;
    moments.mean = origin + mean;
    // The grid spreads the weight of a point uniformly over its step
    moments.variance =
        std::max(sum2 / norm - mean * mean, 0.) + step * step / 12.;

    const double sigma = std::sqrt(moments.variance);
    if (sigma > 4 * step) {
      break;
    }
    low = std::max(low, moments.mean - 10 * sigma);
    high = std::min(high, moments.mean + 10 * sigma);
  }
  return moments;
}

/// The noise-free digital box, optionally constrained by partial edge charges.
/// A single analogue bin uses the digital box: a subthreshold neighbour can
/// remove charge on either side, so its charge need not determine a position.
std::optional<PosteriorMoments> calibrateWithoutNoise(
    const ClusterProjection& p, const ClusterFootprint& footprint,
    double threshold, bool analogue) {
  const double length = footprint.length;
  const double fraction = threshold / footprint.pathLength;
  const double thresholdLength = fraction * length;
  const double minWidth =
      std::min(p.firstUpper - p.firstLower, p.lastUpper - p.lastLower);
  // If a whole neighbour can stay silent, the box would need constraints
  // from further bins. Use the cluster-wide fallback in that case.
  if (fraction >= 1 || (length > 0 && thresholdLength >= minWidth) ||
      (p.beforeWidth > 0 && thresholdLength >= p.beforeWidth) ||
      (p.afterWidth > 0 && thresholdLength >= p.afterWidth)) {
    return std::nullopt;
  }

  // Fired edge bins require an overlap greater than threshold * length /
  // pathLength. Silent neighbours bound the other ends; at a sensor edge,
  // only the hit position must stay inside the sensor.
  const double half = 0.5 * length * (1. - 2. * fraction);
  const double lowFromFirst = p.firstLower + (p.beforeWidth > 0 ? half : 0.);
  const double lowFromLast = p.lastLower - half;
  const double highFromFirst = p.firstUpper + half;
  const double highFromLast = p.lastUpper - (p.afterWidth > 0 ? half : 0.);
  const double low = std::max(lowFromFirst, lowFromLast);
  const double high = std::min(highFromFirst, highFromLast);
  if (low > high) {
    return std::nullopt;
  }

  // The box moves with the footprint length. The centre moves by the mean
  // of the slopes of the active bounds.
  const double halfError = 0.5 * (1. - 2. * fraction) * footprint.lengthError;
  const double centreSlope =
      0.5 *
      ((lowFromFirst >= lowFromLast ? (p.beforeWidth > 0 ? 1. : 0.) : -1.) +
       (highFromFirst <= highFromLast ? 1. : (p.afterWidth > 0 ? -1. : 0.)));
  PosteriorMoments result;
  result.mean = 0.5 * (low + high);
  result.variance =
      (high - low) * (high - low) / 12. + std::pow(centreSlope * halfError, 2);

  // Keep the digital approximation for a single analogue bin.
  if (!analogue || p.firstLower == p.lastLower) {
    return result;
  }

  // A partial edge gives the segment endpoint, hence the track midpoint.
  // A full edge carries no position information. The sign distinguishes the
  // lower endpoint inferred from the first bin from the upper endpoint.
  const double lengthPerCharge = length / footprint.pathLength;
  const double tolerance = 1e-9 * std::max(length, p.lastUpper - p.firstLower);
  const std::array edges = {
      std::array{p.firstUpper, p.firstUpper - p.firstLower, p.firstCharge, 1.},
      std::array{p.lastLower, p.lastUpper - p.lastLower, p.lastCharge, -1.}};
  double sumPosition = 0;
  double sumSlope = 0;
  unsigned int partial = 0;
  for (const auto& [edge, width, charge, sign] : edges) {
    const double overlap = charge * lengthPerCharge;
    if (footprint.lengthError == 0 && overlap > width + tolerance) {
      return std::nullopt;
    }
    if (overlap >= width * (1. - 1e-9)) {
      continue;
    }
    const double position = edge + sign * (0.5 * length - overlap);
    // Contradictory charges must not produce a zero-variance measurement.
    if (footprint.lengthError == 0 &&
        (position < low - tolerance || position > high + tolerance ||
         (partial != 0 && std::abs(position - sumPosition) > tolerance))) {
      return std::nullopt;
    }
    sumPosition += position;
    sumSlope += sign * (0.5 - charge / footprint.pathLength);
    ++partial;
  }
  if (partial != 0) {
    result.mean = std::clamp(sumPosition / partial, low, high);
    result.variance = std::pow(sumSlope / partial * footprint.lengthError, 2);
  }
  return result;
}

/// Calibrate one axis with a uniform position prior and an edge-bin likelihood:
/// first and last bins fire (or have the observed analogue charge), and their
/// touched neighbours stay silent. Uniform deposition and Gaussian cell noise
/// are assumed. With no noise, use the digital box and partial analogue edges;
/// otherwise integrate numerically. If no position is compatible, fall back to
/// a uniform distribution over the cluster.
///
/// Direction uncertainty varies the footprint length at fixed pathLength,
/// ignoring their correlation. This is a small-angle approximation.
PosteriorMoments calibrateClusterAxis(const ClusterProjection& projection,
                                      const ClusterFootprint& footprint,
                                      double threshold, bool analogue,
                                      double chargeNoise) {
  const auto& p = projection;

  // Uniform over the whole cluster if the footprint cannot make the cluster
  const PosteriorMoments fallback{
      0.5 * (p.firstLower + p.lastUpper),
      std::pow(p.lastUpper - p.firstLower, 2) / 12.};

  if (chargeNoise <= 0) {
    return calibrateWithoutNoise(p, footprint, threshold, analogue)
        .value_or(fallback);
  }

  const ClusterLikelihood likelihood{p, footprint.pathLength, threshold,
                                     chargeNoise, analogue};

  const std::array<WeightedLength, 1> central = {{{footprint.length, 1.}}};
  auto moments = integratePosterior(likelihood, p, central);

  // A change of the footprint length moves each end of the segment by half
  // of it. If the Gauss-Hermite nodes move the ends by a part of the width of
  // the posterior, marginalise over the length with a three point rule.
  const double nodeShift = 0.5 * std::sqrt(3.) * footprint.lengthError;
  if (nodeShift > 0 &&
      (!moments || nodeShift * nodeShift > 0.01 * moments->variance)) {
    const std::array<WeightedLength, 3> nodes = {
        {{footprint.length, 2. / 3.},
         {std::max(footprint.length - 2. * nodeShift, 0.), 1. / 6.},
         {footprint.length + 2. * nodeShift, 1. / 6.}}};
    moments = integratePosterior(likelihood, p, nodes);
  }
  return moments.value_or(fallback);
}

}  // namespace

ClusterOnTrackCalibrator::ClusterOnTrackCalibrator(const Config& config)
    : m_cfg(config) {
  if (!std::isfinite(m_cfg.chargeNoise) || m_cfg.chargeNoise < 0 ||
      !std::isfinite(m_cfg.minError) || m_cfg.minError < 0) {
    throw std::invalid_argument(
        "Charge noise and minimum error must be finite and nonnegative");
  }
}

void ClusterOnTrackCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const {
  if (clusters == nullptr) {
    throw std::invalid_argument("ClusterOnTrackCalibrator needs clusters");
  }

  // Start from the measurement as it is
  PassThroughCalibrator{}.calibrate(measurements, clusters, gctx, cctx,
                                    sourceLink, trackState);

  // The fitters allocate the parameters that they calibrate with
  if (!trackState.hasPredicted() && !trackState.hasFiltered() &&
      !trackState.hasSmoothed()) {
    return;
  }
  const auto parameters = trackState.parameters();
  const auto covariance = trackState.covariance();

  const Acts::Surface& surface = trackState.referenceSurface();
  if (surface.type() != Acts::Surface::Plane) {
    return;
  }
  auto digiIt = m_cfg.digitizationConfigs.find(surface.geometryId());
  if (digiIt == m_cfg.digitizationConfigs.end()) {
    return;
  }
  const GeometricConfig& geoCfg = digiIt->geometricDigiConfig;
  if (geoCfg.indices.empty() || geoCfg.segmentation == nullptr ||
      geoCfg.thickness <= 0) {
    return;
  }

  const std::size_t index = sourceLink.get<IndexSourceLink>().index();
  assert(index < clusters->size() && "Cluster index out of bounds");
  const Cluster& cluster = clusters->at(index);
  if (cluster.channels.empty()) {
    return;
  }

  // The direction in the local frame and its derivatives by phi and theta
  const Acts::Vector3 globalDir = Acts::makeDirectionFromPhiTheta(
      parameters[Acts::eBoundPhi], parameters[Acts::eBoundTheta]);
  const Acts::RotationMatrix3 toLocal =
      surface.localToGlobalTransform(gctx).linear().transpose();
  const Acts::Vector3 dir = toLocal * globalDir;
  const Acts::Matrix<3, 2> dirJacobian =
      toLocal * Acts::sphericalToFreeDirectionJacobian(globalDir);
  if (std::abs(dir.z()) < std::numeric_limits<double>::epsilon()) {
    return;
  }
  const Acts::SquareMatrix<2> angleCovariance =
      covariance.block<2, 2>(Acts::eBoundPhi, Acts::eBoundPhi);

  const auto measurement = measurements.getMeasurement(index);
  const std::size_t size = measurement.size();
  Acts::visit_measurement(size, [&](auto N) {
    constexpr std::size_t kSize = decltype(N)::value;
    Acts::Vector<kSize> calibrated = trackState.template calibrated<kSize>();
    Acts::SquareMatrix<kSize> calibratedCovariance =
        trackState.template calibratedCovariance<kSize>();

    for (Acts::BoundIndices bIndex : geoCfg.indices) {
      if ((bIndex != Acts::eBoundLoc0 && bIndex != Acts::eBoundLoc1) ||
          !measurement.contains(bIndex) ||
          static_cast<std::size_t>(bIndex) >= geoCfg.segmentation->getNAxes()) {
        continue;
      }
      const auto axisIndex = static_cast<std::size_t>(bIndex);
      const auto& axis = geoCfg.segmentation->getAxis(bIndex);
      const auto expected = bIndex == Acts::eBoundLoc0
                                ? Acts::AxisDirection::AxisX
                                : Acts::AxisDirection::AxisY;
      if (axis.getDirection() != expected) {
        continue;
      }

      // Project the cluster onto this axis
      std::size_t first = std::numeric_limits<std::size_t>::max();
      std::size_t last = 0;
      for (const auto& cell : cluster.channels) {
        first = std::min<std::size_t>(first, cell.bin[axisIndex]);
        last = std::max<std::size_t>(last, cell.bin[axisIndex]);
      }
      ClusterProjection projection;
      // Cell bins are zero-based while the axis bin indices start at one
      projection.firstLower = axis.getBinLowerBound(first + 1);
      projection.firstUpper = axis.getBinUpperBound(first + 1);
      projection.lastLower = axis.getBinLowerBound(last + 1);
      projection.lastUpper = axis.getBinUpperBound(last + 1);
      if (first > 0) {
        projection.beforeWidth = axis.getBinWidth(first);
      }
      if (last + 1 < axis.getNBins()) {
        projection.afterWidth = axis.getBinWidth(last + 2);
      }
      for (const auto& cell : cluster.channels) {
        if (cell.bin[axisIndex] == first) {
          projection.firstCharge += cell.activation;
          ++projection.firstCells;
        }
        if (cell.bin[axisIndex] == last) {
          projection.lastCharge += cell.activation;
          ++projection.lastCells;
        }
      }

      // The footprint of the predicted track on this axis
      const double slope = dir[bIndex] / dir.z();
      const Acts::Vector2 slopeJacobian =
          (dirJacobian.row(bIndex) - slope * dirJacobian.row(2)).transpose() /
          dir.z();
      ClusterFootprint footprint;
      footprint.length = geoCfg.thickness * std::abs(slope);
      footprint.lengthError =
          geoCfg.thickness *
          std::sqrt(
              std::max(slopeJacobian.dot(angleCovariance * slopeJacobian), 0.));
      footprint.pathLength = geoCfg.thickness / std::abs(dir.z());

      const PosteriorMoments axisCalibration =
          calibrateClusterAxis(projection, footprint, geoCfg.threshold,
                               !geoCfg.digital, m_cfg.chargeNoise);

      const std::size_t i = measurement.indexOf(bIndex);
      calibrated[i] = axisCalibration.mean;
      calibratedCovariance.row(i).setZero();
      calibratedCovariance.col(i).setZero();
      calibratedCovariance(i, i) =
          std::max(axisCalibration.variance, m_cfg.minError * m_cfg.minError);
    }

    trackState.template calibrated<kSize>() = calibrated;
    trackState.template calibratedCovariance<kSize>() = calibratedCovariance;
  });
}

}  // namespace ActsExamples
