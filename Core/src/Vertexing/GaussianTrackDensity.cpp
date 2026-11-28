// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/GaussianTrackDensity.hpp"

#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <span>

namespace Acts {

namespace {
// The bounds already truncate each Gaussian. Index their support, rather than
// testing every track at every trial position. Inserting in input order retains
// the original floating-point summation order in every bin.
class DensityIndex {
 public:
  using Entry = GaussianTrackDensity::TrackEntry;

  explicit DensityIndex(std::span<const Entry> entries) {
    if (entries.size() < 32) {
      return;
    }
    m_lower = std::numeric_limits<double>::infinity();
    m_upper = -std::numeric_limits<double>::infinity();
    for (const auto& entry : entries) {
      if (!std::isfinite(entry.lowerBound) ||
          !std::isfinite(entry.upperBound)) {
        return;
      }
      m_lower = std::min(m_lower, entry.lowerBound);
      m_upper = std::max(m_upper, entry.upperBound);
    }
    const double extent = m_upper - m_lower;
    if (!(extent > 0.) || !std::isfinite(extent)) {
      return;
    }
    const std::size_t nBins = entries.size();
    m_scale = static_cast<double>(nBins) / extent;
    if (!std::isfinite(m_scale)) {
      return;
    }
    m_offsets.resize(nBins + 1);
    // Wide, overlapping supports can make an index counterproductive. Bound
    // its memory/work and use the original scan in that case.
    std::size_t remaining = 32 * entries.size();
    for (const auto& entry : entries) {
      if (!(entry.lowerBound < entry.upperBound)) {
        continue;
      }
      const auto first = bin(entry.lowerBound);
      const auto last = bin(entry.upperBound);
      const auto count = last - first + 1;
      if (count > remaining) {
        m_offsets.clear();
        return;
      }
      remaining -= count;
      for (std::size_t i = first; i <= last; ++i) {
        ++m_offsets[i + 1];
      }
    }
    std::partial_sum(m_offsets.begin(), m_offsets.end(), m_offsets.begin());
    m_entries.resize(m_offsets.back());
    for (const auto& entry : entries) {
      if (!(entry.lowerBound < entry.upperBound)) {
        continue;
      }
      const auto first = bin(entry.lowerBound);
      const auto last = bin(entry.upperBound);
      for (std::size_t i = first; i <= last; ++i) {
        m_entries[m_offsets[i]++] = &entry;
      }
    }
    // The insertion cursors now hold the end of each bin. Restore its start.
    for (std::size_t i = nBins; i > 0; --i) {
      m_offsets[i] = m_offsets[i - 1];
    }
    m_offsets[0] = 0;
  }

  // nullopt selects the original scan; an empty span means no finite support.
  std::optional<std::span<const Entry* const>> candidates(double z) const {
    if (m_offsets.empty()) {
      return std::nullopt;
    }
    if (!std::isfinite(z) || z < m_lower || z > m_upper) {
      return std::span<const Entry* const>{};
    }
    const auto i = bin(z);
    return std::span<const Entry* const>{m_entries}.subspan(
        m_offsets[i], m_offsets[i + 1] - m_offsets[i]);
  }

 private:
  std::size_t bin(double z) const {
    return std::min(static_cast<std::size_t>((z - m_lower) * m_scale),
                    m_offsets.size() - 2);
  }

  double m_lower = 0.;
  double m_upper = 0.;
  double m_scale = 0.;
  std::vector<std::size_t> m_offsets;
  std::vector<const Entry*> m_entries;
};
}  // namespace

Result<std::optional<std::pair<double, double>>>
Acts::GaussianTrackDensity::globalMaximumWithWidth(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto result = addTracks(state, trackList);
  if (!result.ok()) {
    return result.error();
  }

  const DensityIndex index(state.trackEntries);
  const auto densityAt = [&](double z) {
    const auto candidates = index.candidates(z);
    if (!candidates) {
      return trackDensityAndDerivatives(state, z);
    }
    GaussianTrackDensityStore density(z);
    for (const auto* entry : *candidates) {
      density.addTrackToDensity(*entry);
    }
    return density.densityAndDerivatives();
  };

  double maxPosition = 0.;
  double maxDensity = 0.;
  double maxSecondDerivative = 0.;

  for (const auto& track : state.trackEntries) {
    double trialZ = track.z;

    auto [density, firstDerivative, secondDerivative] = densityAt(trialZ);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);

    trialZ += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) = densityAt(trialZ);

    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) = densityAt(trialZ);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
  }

  if (maxSecondDerivative == 0.) {
    return std::nullopt;
  }

  return std::pair{maxPosition, std::sqrt(-(maxDensity / maxSecondDerivative))};
}

Result<std::optional<double>> Acts::GaussianTrackDensity::globalMaximum(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto maxRes = globalMaximumWithWidth(state, trackList);
  if (!maxRes.ok()) {
    return maxRes.error();
  }
  const auto& maxOpt = *maxRes;
  if (!maxOpt.has_value()) {
    return std::nullopt;
  }
  return maxOpt->first;
}

Result<void> Acts::GaussianTrackDensity::addTracks(
    State& state, const std::vector<InputTrack>& trackList) const {
  for (auto trk : trackList) {
    const BoundTrackParameters& boundParams = m_cfg.extractParameters(trk);
    // Get required track parameters
    const double d0 = boundParams.parameters()[BoundIndices::eBoundLoc0];
    const double z0 = boundParams.parameters()[BoundIndices::eBoundLoc1];
    // Get track covariance
    if (!boundParams.covariance().has_value()) {
      return VertexingError::NoCovariance;
    }
    const auto perigeeCov = *(boundParams.covariance());
    const double covDD =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0);
    const double covZZ =
        perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1);
    const double covDZ =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1);
    const double covDeterminant = (perigeeCov.block<2, 2>(0, 0)).determinant();

    // Do track selection based on track cov matrix and m_cfg.d0SignificanceCut
    if ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) ||
        (covZZ <= 0) || (covDeterminant <= 0)) {
      continue;
    }

    // Calculate track density quantities
    double constantTerm =
        -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ) /
        (2. * covDeterminant);
    const double linearTerm =
        (d0 * covDZ + z0 * covDD) /
        covDeterminant;  // minus signs and factors of 2 cancel...
    const double quadraticTerm = -covDD / (2. * covDeterminant);
    double discriminant =
        linearTerm * linearTerm -
        4. * quadraticTerm * (constantTerm + 2. * m_cfg.z0SignificanceCut);
    if (discriminant < 0) {
      continue;
    }

    // Add the track to the current maps in the state
    discriminant = std::sqrt(discriminant);
    const double zMax = (-linearTerm - discriminant) / (2. * quadraticTerm);
    const double zMin = (-linearTerm + discriminant) / (2. * quadraticTerm);
    constantTerm -= std::log(2. * std::numbers::pi * std::sqrt(covDeterminant));

    state.trackEntries.emplace_back(z0, constantTerm, linearTerm, quadraticTerm,
                                    zMin, zMax);
  }
  return Result<void>::success();
}

std::tuple<double, double, double>
Acts::GaussianTrackDensity::trackDensityAndDerivatives(State& state,
                                                       double z) const {
  GaussianTrackDensityStore densityResult(z);
  for (const auto& trackEntry : state.trackEntries) {
    densityResult.addTrackToDensity(trackEntry);
  }
  return densityResult.densityAndDerivatives();
}

std::tuple<double, double, double> Acts::GaussianTrackDensity::updateMaximum(
    double newZ, double newValue, double newSecondDerivative, double maxZ,
    double maxValue, double maxSecondDerivative) const {
  if (newValue > maxValue) {
    maxZ = newZ;
    maxValue = newValue;
    maxSecondDerivative = newSecondDerivative;
  }
  return {maxZ, maxValue, maxSecondDerivative};
}

double Acts::GaussianTrackDensity::stepSize(double y, double dy,
                                            double ddy) const {
  return (m_cfg.isGaussianShaped ? (y * dy) / (dy * dy - y * ddy) : -dy / ddy);
}

void Acts::GaussianTrackDensity::GaussianTrackDensityStore::addTrackToDensity(
    const TrackEntry& entry) {
  // Take track only if it's within bounds
  if (entry.lowerBound < m_z && m_z < entry.upperBound) {
    double delta = std::exp(entry.c0 + m_z * (entry.c1 + m_z * entry.c2));
    double qPrime = entry.c1 + 2. * m_z * entry.c2;
    double deltaPrime = delta * qPrime;
    m_density += delta;
    m_firstDerivative += deltaPrime;
    m_secondDerivative += 2. * entry.c2 * delta + qPrime * deltaPrime;
  }
}

}  // namespace Acts
