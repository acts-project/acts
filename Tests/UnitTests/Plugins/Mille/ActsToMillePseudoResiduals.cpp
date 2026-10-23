// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Mille/MilleDecoder.h"
#include "Mille/MilleFactory.h"

// The Kalman correlations are written as pseudo-measurements along the
// eigenvectors of the correlation term. The record is expanded around the
// smoothed track, whose measurement chi2 gradient H^T V^-1 r is in general not
// zero (multiple scattering kinks). The pseudo-measurements must carry
// residuals that cancel it, so that the local fit in Millepede is at its
// minimum on the smoothed track. Otherwise the local fit moves the track and
// the alignment derivatives are biased.

BOOST_AUTO_TEST_SUITE(ActsToMillePseudoResidualsTests)

BOOST_AUTO_TEST_CASE(SmoothedTrackIsLocalFitMinimum) {
  // three track states, loc0 and loc1 measured on each
  constexpr std::size_t nStates = 3;
  constexpr std::size_t nPar = nStates * Acts::eBoundSize;
  constexpr std::size_t nMeas = 2 * nStates;
  // the last measurement is on a surface that is not aligned
  constexpr std::size_t iMeasNotAligned = nMeas - 1;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::RectangleBounds>(10., 10.));

  ActsAlignment::detail::TrackAlignmentState state;
  state.measurementDim = nMeas;
  state.trackParametersDim = nPar;
  state.alignmentDof = Acts::eAlignmentSize;
  state.alignedSurfaces.emplace(surface.get(), std::make_pair(0u, 0u));

  state.measurementCovariance =
      Acts::DynamicMatrix::Identity(nMeas, nMeas) * 0.05 * 0.05;
  state.projectionMatrix = Acts::DynamicMatrix::Zero(nMeas, nPar);
  for (std::size_t s = 0; s < nStates; ++s) {
    const std::size_t o = s * Acts::eBoundSize;
    state.projectionMatrix(2 * s, o + Acts::eBoundLoc0) = 1.;
    state.projectionMatrix(2 * s + 1, o + Acts::eBoundLoc1) = 1.;
  }

  // Track parameter covariance of a Kalman fit: the inverse of the measurement
  // information plus a dense, positive definite correlation term K, which
  // couples the parameters on different surfaces.
  Acts::DynamicMatrix b(nPar, nPar);
  for (std::size_t i = 0; i < nPar; ++i) {
    for (std::size_t j = 0; j < nPar; ++j) {
      b(i, j) = std::sin(1. + i + 2. * j);
    }
  }
  const Acts::DynamicMatrix correlationTerm =
      b.transpose() * b + Acts::DynamicMatrix::Identity(nPar, nPar);
  const Acts::DynamicMatrix weightMatMeasurements =
      state.projectionMatrix.transpose() *
      state.measurementCovariance.inverse() * state.projectionMatrix;
  state.trackParametersCovariance =
      (weightMatMeasurements + correlationTerm).inverse();

  // residuals of a track with kinks: the measurement chi2 gradient is not zero
  state.residual = Acts::DynamicVector::Zero(nMeas);
  state.alignmentToResidualDerivative =
      Acts::DynamicMatrix::Zero(nMeas, state.alignmentDof);
  for (std::size_t i = 0; i < nMeas; ++i) {
    state.residual(i) = 0.01 * (i % 2 == 0 ? 1. : -2.) * (i + 1);
    if (i == iMeasNotAligned) {
      continue;
    }
    for (std::size_t a = 0; a < Acts::eAlignmentSize; ++a) {
      state.alignmentToResidualDerivative(i, a) = 0.1 * (a + 1) + 0.01 * i;
    }
  }
  const Acts::DynamicVector measGradient =
      state.projectionMatrix.transpose() *
      state.measurementCovariance.inverse() * state.residual;
  BOOST_REQUIRE(measGradient.norm() > 1.);

  auto logger = Acts::getDefaultLogger("ActsToMillePseudoResiduals",
                                       Acts::Logging::WARNING);

  const std::string fname = "ActsToMillePseudoResiduals.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    ActsPlugins::ActsToMille::dumpToMille(state, *out, false, *logger);
  }

  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));
  Mille::MilleDecoder decoder;
  std::vector<Mille::MilleMeasurement> measurements;
  BOOST_REQUIRE(decoder.decode(*reader, measurements) ==
                Mille::MilleDecoder::ReadResult::OK);
  BOOST_REQUIRE_EQUAL(measurements.size(), nMeas + nPar);

  // Gradient of the local fit chi2 at the expansion point (the smoothed
  // track), summed over the measurements and the pseudo-measurements. Each
  // term is compared to the scale of the measurement contribution; the record
  // is written in single precision.
  Acts::DynamicVector localGradient = Acts::DynamicVector::Zero(nPar);
  std::size_t nPseudoNonZero = 0;
  for (const auto& measurement : measurements) {
    const double weight =
        1. / (measurement.uncertainty * measurement.uncertainty);
    for (std::size_t iLoc = 0; iLoc < measurement.localLabels.size(); ++iLoc) {
      localGradient(measurement.localLabels[iLoc] - 1) +=
          weight * measurement.measurement * measurement.localDerivatives[iLoc];
    }
    if (measurement.globalLabels.empty() &&
        measurement.localLabels.size() > 1 && measurement.measurement != 0) {
      ++nPseudoNonZero;
    }
  }
  BOOST_CHECK_EQUAL(nPseudoNonZero, nPar);
  BOOST_CHECK_SMALL(localGradient.norm() / measGradient.norm(), 1e-5);

  // The reader still separates the surface measurements from the
  // pseudo-measurements, including the one without alignment derivatives.
  auto reader2 = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader2 != nullptr);
  BOOST_REQUIRE(reader2->open(fname));
  const std::unordered_map<const Acts::Surface*, std::size_t>
      idxedAlignSurfaces{{surface.get(), 0}};
  ActsAlignment::detail::TrackAlignmentState readState;
  BOOST_REQUIRE(ActsPlugins::ActsToMille::unpackMilleRecord(
                    *reader2, readState, idxedAlignSurfaces, *logger) ==
                Mille::MilleDecoder::ReadResult::OK);
  BOOST_CHECK_EQUAL(readState.measurementDim, nMeas);
  BOOST_CHECK_EQUAL(readState.trackParametersDim, nPar);
  BOOST_REQUIRE_EQUAL(readState.residual.size(), nMeas);
  for (std::size_t i = 0; i < nMeas; ++i) {
    BOOST_CHECK_CLOSE(readState.residual(i), state.residual(i), 1e-4);
  }
  BOOST_CHECK(readState.projectionMatrix.isApprox(state.projectionMatrix));
}

BOOST_AUTO_TEST_SUITE_END()
