// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/GaussianTrackDensity.hpp"

#include "Acts/Vertexing/VertexingError.hpp"

#include <math.h>

namespace Acts {

Result<std::optional<std::pair<double, double>>>
Acts::GaussianTrackDensity::globalMaximumWithWidth(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto result = addTracks(state, trackList);
  if (!result.ok()) {
    return result.error();
  }

  double maxPosition = 0.;
  double maxDensity = 0.;
  double maxSecondDerivative = 0.;

//stand in variables
  double other1=0.;
  double other2=0.;
  double other3=0.;
  double other4=0.;
  double other5=0.;
  double other6=0.;

if (state.trackEntries.empty()) {
    std::cout << "Track Entries is empty." << std::endl;}
  for (const auto& track : state.trackEntries) {
    double trialZ = track.z;
    double time = track.time;

    //,tfirstDerivative, tsecondDerivative,dzdt] =

    auto [density, zfirstDerivative, zsecondDerivative,tfirstDerivative, tsecondDerivative,dzdt] =
        trackDensityAndDerivatives(state, trialZ,time);

    std::cout << "Truth statement: " << (zsecondDerivative >= 0. || density <= 0.) << std::endl;
    if (zsecondDerivative >= 0. || density <= 0.) {
      continue;
    }
    //added other1-3 to have 6 params
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, zsecondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);

    trialZ += stepSize(density, zfirstDerivative, zsecondDerivative);
    std::tie(density, zfirstDerivative, zsecondDerivative, other1, other2, other3) =
        trackDensityAndDerivatives(state, trialZ,time);

    if (zsecondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, zsecondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
    trialZ += stepSize(density, zfirstDerivative, zsecondDerivative);
    std::tie(density, zfirstDerivative, zsecondDerivative,other4,other5,other6) =
        trackDensityAndDerivatives(state, trialZ,time);
    if (zsecondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, zsecondDerivative, maxPosition,
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
    const double time = boundParams.parameters()[BoundIndices::eBoundTime];
    //const double time = InputTrack::extractParameters(trackList[0]).time();

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

    //Add time into Cov 
    const double covDT = perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundTime); // = 0
    const double covZT = perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundTime); // = 0
    const double covTT = perigeeCov(BoundIndices::eBoundTime, BoundIndices::eBoundTime); // = 1 ?


   //const double covDeterminant = (perigeeCov.block<2, 2>(0, 0)).determinant();

    // Update determinant for 3x3 matrix because I'm not sure how perigee works
    Eigen::Matrix3d covMatrix;
    covMatrix << covDD, covDZ, covDT,
                 covDZ, covZZ, covZT,
                 covDT, covZT, covTT;
    const double covDeterminant = covMatrix.determinant();

    //std::cout << "covDD: " << covDD << std::endl;
    //std::cout << "covZT: " << covZZ << std::endl;
    //std::cout << "covDeterminant: " << covDeterminant << std::endl;

    // Do track selection based on track cov matrix and m_cfg.d0SignificanceCut
    if ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) ||
        (covZZ <= 0) || (covDeterminant <= 0)) {
          std::cout << "Truth statement: " << ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) || (covZZ <= 0) || (covDeterminant <= 0)) << std::endl;
      continue;
    }

    //std::cout << "1st part of truth: " << ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut)) << std::endl;
    //std::cout << "2nd part of truth: " << ((covZZ <= 0) || (covDeterminant <= 0)) << std::endl;
    //std::cout << "Truth statement: " << ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) ||
        //(covZZ <= 0) || (covDeterminant <= 0)) << std::endl;



    // Calculate track density quantities
    //Change time to covTT
    double constantTerm =
        -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ +
          2. * d0 * time * covDT + 2. * z0 * time * covZT + time * time * covTT) /
        (2. * covDeterminant); //has to be 0 or less

    // double constantTerm =
    //     -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ) /
    //     (2. * covDeterminant);
    const double linearTerm =
        (d0 * covDZ + z0 * covDD + time * covZT) / covDeterminant;
    const double quadraticTerm = -covDD / (2. * covDeterminant);
    double discriminant =(
        linearTerm * linearTerm -
        4. * quadraticTerm * (constantTerm + 2. * m_cfg.z0SignificanceCut));

    std::cout << "Constant Term: " << constantTerm << std::endl;
    std::cout << "d0: " << d0 << std::endl;
    std::cout << "z0: " << z0 << std::endl;
    std::cout << "time: " << time << std::endl;
    std::cout << "covZZ: " << covZZ << std::endl;
    std::cout << "covDD: " << covDD << std::endl;
    std::cout << "covDZ: " << covDZ << std::endl;
    std::cout << "covDT: " << covDT << std::endl;
    std::cout << "covZT: " << covZT << std::endl;
    std::cout << "covTT: " << covTT << std::endl;
    std::cout << "covDeterminant: " << covDeterminant << std::endl;
    std::cout << "discriminant: " << (discriminant) << std::endl;  

    if (discriminant < 0) {
      
      continue; //continue means skip the rest i think --> discrimant should be higher than 0!
    }
  
    // Add the track to the current maps in the state
    discriminant = std::sqrt(discriminant);
    const double zMax = (-linearTerm - discriminant) / (2. * quadraticTerm);
    const double zMin = (-linearTerm + discriminant) / (2. * quadraticTerm);
    constantTerm -= std::log(2. * M_PI * std::sqrt(covDeterminant));

    state.trackEntries.emplace_back(z0, constantTerm, linearTerm, quadraticTerm,
                                    zMin, zMax,time); // Included time

  }
  return Result<void>::success();
}

std::tuple<double, double, double, double, double, double> //added 3 doubles since output has 6 'param' --> solved something (not what I was aiming for but less errors)
Acts::GaussianTrackDensity::trackDensityAndDerivatives(State& state,
                                                       double z, double time) const {
  GaussianTrackDensityStore densityResult(z,time,m_cfg.t0SignificanceCut);
  for (const auto& trackEntry : state.trackEntries) {
    std::cout << "Stuff: " << std::endl;  
    densityResult.addTrackToDensity(trackEntry);
  }
  auto [density, zfirstDerivative, zsecondDerivative, tfirstDerivative, tsecondDerivative, dzdt] = densityResult.densityAndDerivatives();

  std::cout << "Density: " << density 
          << ", zFirstDerivative: " << zfirstDerivative 
          << ", zSecondDerivative: " << zsecondDerivative 
          << ", tFirstDerivative: " << tfirstDerivative 
          << ", tSecondDerivative: " << tsecondDerivative 
          << ", dz/dt: " << dzdt 
          << std::endl;

  return densityResult.densityAndDerivatives(); //what does this even return?
}

std::tuple<double, double, double> Acts::GaussianTrackDensity::updateMaximum( //added 3 doubles and tests to return same size-ish --> more errors and not helping --> if max time then need to implement
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
    m_zfirstDerivative += deltaPrime;
    m_zsecondDerivative += 2. * entry.c2 * delta + qPrime * deltaPrime;
    m_tfirstDerivative += delta;
    m_tsecondDerivative += delta;
    m_dzdt += delta;
  
  }
}

// double Acts::GaussianTrackDensity::GaussianTrackDensityStore::computeTimeFactor(
//     double entryTime, double currentTime) const {
//   return std::exp(-std::abs(entryTime - currentTime) / m_cfg.timeScale);

// } 
} // namespace Acts