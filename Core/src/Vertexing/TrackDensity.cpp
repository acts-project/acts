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
  // Get required track parameters
  const double d0 = trk.parameters()[ParID_t::eLOC_D0];
  const double z0 = trk.parameters()[ParID_t::eLOC_Z0];
  // Get track covariance
  const auto perigeeCov = *(trk.covariance());
  const double covDD = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0);
  const double covZZ = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eLOC_Z0);
  const double covDZ = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_Z0);
  const double covDeterminant = covDD * covZZ - covDZ * covDZ;

  // Do track selection based on track cov matrix and d0SignificanceCut
  if ((covDD <= 0) || (d0 * d0 / covDD > d0SignificanceCut) || (covZZ <= 0) ||
      (covDeterminant <= 0)) {
    return;
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
      4. * quadraticTerm * (constantTerm + 2. * z0SignificanceCut);
  if (discriminant < 0) {
    return;
  }

  // Add the track to the current maps in the state
  discriminant = std::sqrt(discriminant);
  const double zMax = (-linearTerm - discriminant) / (2. * quadraticTerm);
  const double zMin = (-linearTerm + discriminant) / (2. * quadraticTerm);
  constantTerm -= std::log(2. * M_PI * std::sqrt(covDeterminant));
  state.trackEntries.emplace_back(z0, constantTerm, linearTerm, quadraticTerm,
                                  zMin, zMax);
}

std::pair<double, double> Acts::TrackDensity::globalMaximumWithWidth(
    State& state) const {
  double maxPosition = 0.;
  double maxDensity = 0.;
  double maxSecondDerivative = 0.;

  for (const auto& track : state.trackEntries) {
    double trialZ = track.z;
    double density = 0.;
    double firstDerivative = 0.;
    double secondDerivative = 0.;
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
  }

  return std::make_pair(maxPosition,
                        std::sqrt(-(maxDensity / maxSecondDerivative)));
}

double Acts::TrackDensity::globalMaximum(State& state) const {
  return globalMaximumWithWidth(state).first;
}

double Acts::TrackDensity::trackDensity(State& state, double z) const {
  double firstDerivative = 0;
  double secondDerivative = 0;
  return trackDensity(state, z, firstDerivative, secondDerivative);
}

double Acts::TrackDensity::trackDensity(State& state, double z,
                                        double& firstDerivative,
                                        double& secondDerivative) const {
  TrackDensityStore densityResult(z);
  for (const auto& trackEntry : state.trackEntries) {
    densityResult.addTrackToDensity(trackEntry);
  }
  firstDerivative = densityResult.firstDerivative();
  secondDerivative = densityResult.secondDerivative();

  return densityResult.density();
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

void Acts::TrackDensity::TrackDensityStore::addTrackToDensity(
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
