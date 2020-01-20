// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/TrackDensity.hpp"
#include <math.h>

void Acts::TrackDensity::addTrack(State& state, const BoundParameters& trk,
                                  const double d0SignificanceCut,
                                  const double z0SignificanceCut) const {
  if (state.trackMap.count(trk) != 0) {
    return;
  }
  const double d0 = trk.parameters()[ParID_t::eLOC_D0];
  const double z0 = trk.parameters()[ParID_t::eLOC_Z0];
  const auto perigeeCov = *(trk.covariance());
  const double covDD = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0);
  const double covZZ = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eLOC_Z0);

  const double covDZ = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_Z0);

  const double covDeterminant = covDD * covZZ - covDZ * covDZ;

  if ((covDD <= 0) || (d0 * d0 / covDD > d0SignificanceCut) || (covZZ <= 0) ||
      (covDeterminant <= 0)) {
    return;
  }

  double constantTerm =
      -(d0 * d0 * covZZ + z0 * z0 * covDD + 2 * d0 * z0 * covDZ) /
      (2 * covDeterminant);
  const double linearTerm =
      (d0 * covDZ + z0 * covDD) /
      covDeterminant;  // minus signs and factors of 2 cancel...
  const double quadraticTerm = -covDD / (2 * covDeterminant);
  double discriminant =
      linearTerm * linearTerm -
      4 * quadraticTerm * (constantTerm + 2 * z0SignificanceCut);
  if (discriminant < 0) {
    return;
  }

  discriminant = std::sqrt(discriminant);
  const double zMax = (-linearTerm - discriminant) / (2 * quadraticTerm);
  const double zMin = (-linearTerm + discriminant) / (2 * quadraticTerm);
  state.maxZRange = std::max(state.maxZRange, std::max(zMax - z0, z0 - zMin));
  constantTerm -= std::log(2 * M_PI * std::sqrt(covDeterminant));
  state.trackMap.emplace(std::piecewise_construct, std::forward_as_tuple(trk),
                         std::forward_as_tuple(constantTerm, linearTerm,
                                               quadraticTerm, zMin, zMax));
  state.lowerMap.emplace(std::piecewise_construct,
                         std::forward_as_tuple(constantTerm, linearTerm,
                                               quadraticTerm, zMin, zMax),
                         std::forward_as_tuple(trk));
  state.upperMap.emplace(std::piecewise_construct,
                         std::forward_as_tuple(constantTerm, linearTerm,
                                               quadraticTerm, zMin, zMax),
                         std::forward_as_tuple(trk));
}

std::pair<double, double> Acts::TrackDensity::globalMaximumWithWidth(
    State& state) const {
  double maximumPosition = 0.0;
  double maximumDensity = 0.0;
  double maxCurvature = 0.;

  for (const auto& entry : state.trackMap) {
    double trialZ = entry.first.parameters()[ParID_t::eLOC_Z0];
    double density = 0.0;
    double slope = 0.0;
    double curvature = 0.0;
    density = trackDensity(state, trialZ, slope, curvature);
    if (curvature >= 0.0 || density <= 0.0) {
      continue;
    }
    updateMaximum(trialZ, density, curvature, maximumPosition, maximumDensity,
                  maxCurvature);
    trialZ += stepSize(density, slope, curvature);
    density = trackDensity(state, trialZ, slope, curvature);
    if (curvature >= 0.0 || density <= 0.0) {
      continue;
    }
    updateMaximum(trialZ, density, curvature, maximumPosition, maximumDensity,
                  maxCurvature);
    trialZ += stepSize(density, slope, curvature);
    density = trackDensity(state, trialZ, slope, curvature);
    if (curvature >= 0.0 || density <= 0.0) {
      continue;
    }
    updateMaximum(trialZ, density, curvature, maximumPosition, maximumDensity,
                  maxCurvature);
  }
  if (maximumDensity <= 0) {
    std::cout << "Global maximum at density of 0; track map contains "
              << state.trackMap.size() << " tracks" << std::endl;
  }

  return std::make_pair(maximumPosition,
                        std::sqrt(-(maximumDensity / maxCurvature)));
}

double Acts::TrackDensity::globalMaximum(State& state) const {
  return globalMaximumWithWidth(state).first;
}

double Acts::TrackDensity::trackDensity(State& state, double z) const {
  double sum = 0.0;

  TrackEntry target(z);
  TrackMap overlaps;
  LowerMap::const_iterator left = state.lowerMap.lower_bound(
      target);  // first track whose UPPER bound is not less than z
  if (left == state.lowerMap.end()) {
    // z is to the right of every track's range
    return sum;
  }

  UpperMap::const_iterator right = state.upperMap.upper_bound(
      target);  // first track whose LOWER bound is greater than z
  if (right == state.upperMap.begin()) {
    // z is to the left of every track's range
    return sum;
  }

  for (auto itrk = left; itrk != state.lowerMap.end(); itrk++) {
    if (itrk->first.upperBound > z + state.maxZRange) {
      break;
    }
    if (z >= itrk->first.lowerBound && z <= itrk->first.upperBound) {
      overlaps[itrk->second] = itrk->first;
    }
  }
  for (auto itrk = right; itrk-- != state.upperMap.begin();) {
    if (itrk->first.lowerBound < z - state.maxZRange) {
      break;
    }
    if (z >= itrk->first.lowerBound && z <= itrk->first.upperBound) {
      overlaps[itrk->second] = itrk->first;
    }
  }
  for (const auto& entry : overlaps) {
    sum +=
        std::exp(entry.second.c0 + z * (entry.second.c1 + z * entry.second.c2));
  }
  return sum;
}

double Acts::TrackDensity::trackDensity(State& state, double z,
                                        double& firstDerivative,
                                        double& secondDerivative) const {
  double density = 0.0;
  firstDerivative = 0.0;
  secondDerivative = 0.0;
  TrackEntry target(z);
  TrackMap overlaps;

  LowerMap::const_iterator left = state.lowerMap.lower_bound(
      target);  // first track whose UPPER bound is not less than z
  if (left == state.lowerMap.end()) {
    // z is to the right of every track's range
    return density;
  }

  UpperMap::const_iterator right = state.upperMap.upper_bound(
      target);  // first track whose LOWER bound is greater than z
  if (right == state.upperMap.begin()) {
    // z is to the left of every track's range
    return density;
  }

  for (auto itrk = left; itrk != state.lowerMap.end(); itrk++) {
    if (itrk->first.upperBound > z + state.maxZRange) {
      break;
    }
    if (z >= itrk->first.lowerBound && z <= itrk->first.upperBound) {
      overlaps[itrk->second] = itrk->first;
    }
  }
  for (auto itrk = right; itrk-- != state.upperMap.begin();) {
    if (itrk->first.lowerBound < z - state.maxZRange) {
      break;
    }
    if (z >= itrk->first.lowerBound && z <= itrk->first.upperBound) {
      overlaps[itrk->second] = itrk->first;
    }
  }

  for (const auto& entry : overlaps) {
    if (entry.second.lowerBound > z || entry.second.upperBound < z) {
      continue;
    }
    double delta =
        std::exp(entry.second.c0 + z * (entry.second.c1 + z * entry.second.c2));
    density += delta;
    double qPrime = entry.second.c1 + 2 * z * entry.second.c2;
    double deltaPrime = delta * qPrime;
    firstDerivative += deltaPrime;
    secondDerivative += 2 * entry.second.c2 * delta + qPrime * deltaPrime;
  }

  return density;
}

void Acts::TrackDensity::updateMaximum(double newZ, double newValue,
                                       double newSecondDerivative, double& maxZ,
                                       double& maxValue,
                                       double& maxSecondDerivative) const {
  if (newValue > maxValue) {
    maxZ = newZ;
    maxValue = newValue;
    maxSecondDerivative = newSecondDerivative;
  }
}

double Acts::TrackDensity::stepSize(double y, double dy, double ddy) const {
  return (m_cfg.isGaussianShaped ? (y * dy) / (dy * dy - y * ddy) : -dy / ddy);
}
