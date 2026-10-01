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

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Mille/MilleDecoder.h"
#include "Mille/MilleFactory.h"

// removeUnconstrainedTrackPar drops track parameters whose variance jumps by
// more than 1e6 over the smaller ones. The variances have different units, so
// the measured positions can sit above such a jump. They must be kept: a
// measurement written without its local derivatives would constrain the
// alignment as if the track could not move.

BOOST_AUTO_TEST_SUITE(ActsToMilleDiscardGuardTests)

BOOST_AUTO_TEST_CASE(MeasuredParametersAreKept) {
  // two track states, loc0 and loc1 measured on each
  constexpr std::size_t nStates = 2;
  constexpr std::size_t nPar = nStates * Acts::eBoundSize;
  constexpr std::size_t nMeas = 2 * nStates;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::RectangleBounds>(10., 10.));

  ActsAlignment::detail::TrackAlignmentState state;
  state.measurementDim = nMeas;
  state.trackParametersDim = nPar;
  state.alignmentDof = Acts::eAlignmentSize;
  state.alignedSurfaces.emplace(surface.get(), std::make_pair(0u, 0u));

  // measurement sigma 50 um; the smoothed positions are better known than
  // that (variance 1e-3 mm^2), the angles much better (1e-12 rad^2): the
  // positions sit 1e9 above the angles and the old criterion dropped them.
  // The time is unconstrained and is dropped in any case.
  state.measurementCovariance =
      Acts::DynamicMatrix::Identity(nMeas, nMeas) * 0.05 * 0.05;
  state.trackParametersCovariance = Acts::DynamicMatrix::Zero(nPar, nPar);
  state.projectionMatrix = Acts::DynamicMatrix::Zero(nMeas, nPar);
  for (std::size_t s = 0; s < nStates; ++s) {
    const std::size_t o = s * Acts::eBoundSize;
    state.trackParametersCovariance(o + Acts::eBoundLoc0,
                                    o + Acts::eBoundLoc0) = 1e-3;
    state.trackParametersCovariance(o + Acts::eBoundLoc1,
                                    o + Acts::eBoundLoc1) = 1e-3;
    state.trackParametersCovariance(o + Acts::eBoundPhi, o + Acts::eBoundPhi) =
        1e-12;
    state.trackParametersCovariance(o + Acts::eBoundTheta,
                                    o + Acts::eBoundTheta) = 1e-12;
    state.trackParametersCovariance(o + Acts::eBoundQOverP,
                                    o + Acts::eBoundQOverP) = 1e-10;
    state.trackParametersCovariance(o + Acts::eBoundTime,
                                    o + Acts::eBoundTime) = 1e2;
    state.projectionMatrix(2 * s, o + Acts::eBoundLoc0) = 1.;
    state.projectionMatrix(2 * s + 1, o + Acts::eBoundLoc1) = 1.;
  }
  state.residual = Acts::DynamicVector::Zero(nMeas);
  state.alignmentToResidualDerivative =
      Acts::DynamicMatrix::Zero(nMeas, state.alignmentDof);
  for (std::size_t i = 0; i < nMeas; ++i) {
    state.residual(i) = 0.01 * (i + 1);
    for (std::size_t a = 0; a < Acts::eAlignmentSize; ++a) {
      state.alignmentToResidualDerivative(i, a) = 0.1 * (a + 1) + 0.01 * i;
    }
  }

  std::ostringstream logStream;
  auto logger = Acts::getDefaultLogger("ActsToMilleDiscardGuard",
                                       Acts::Logging::WARNING, &logStream);

  const std::string fname = "ActsToMilleDiscardGuard.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    ActsPlugins::ActsToMille::dumpToMille(state, *out, true, *logger);
  }
  BOOST_CHECK(logStream.str().find("measured track parameter") !=
              std::string::npos);

  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));
  Mille::MilleDecoder decoder;
  std::vector<Mille::MilleMeasurement> measurements;
  BOOST_REQUIRE(decoder.decode(*reader, measurements) ==
                Mille::MilleDecoder::ReadResult::OK);

  // every measurement with alignment derivatives keeps a local derivative
  std::size_t nWithGlobals = 0;
  for (const auto& measurement : measurements) {
    if (measurement.globalLabels.empty()) {
      continue;
    }
    ++nWithGlobals;
    BOOST_CHECK(!measurement.localLabels.empty());
  }
  BOOST_CHECK_EQUAL(nWithGlobals, nMeas);
}

BOOST_AUTO_TEST_SUITE_END()
