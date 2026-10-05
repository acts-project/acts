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
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "Mille/MilleDecoder.h"
#include "Mille/MilleFactory.h"

// dumpToMille is called from several threads into the same Mille record, e.g.
// by an algorithm running on several events at once. Every track must end up
// in a record of its own: two tracks merged into one record are fitted as a
// single track by Millepede.

namespace {

constexpr std::size_t nStates = 2;
constexpr std::size_t nPar = nStates * Acts::eBoundSize;
constexpr std::size_t nMeas = 2 * nStates;

/// A track state whose first residual is its track number, so that every
/// record can be traced back to the track it was written from.
ActsAlignment::detail::TrackAlignmentState makeState(
    const Acts::Surface& surface, std::size_t track) {
  ActsAlignment::detail::TrackAlignmentState state;
  state.measurementDim = nMeas;
  state.trackParametersDim = nPar;
  state.alignmentDof = Acts::eAlignmentSize;
  state.alignedSurfaces.emplace(&surface, std::make_pair(0u, 0u));

  state.measurementCovariance =
      Acts::DynamicMatrix::Identity(nMeas, nMeas) * 0.05 * 0.05;
  state.projectionMatrix = Acts::DynamicMatrix::Zero(nMeas, nPar);
  for (std::size_t s = 0; s < nStates; ++s) {
    const std::size_t o = s * Acts::eBoundSize;
    state.projectionMatrix(2 * s, o + Acts::eBoundLoc0) = 1.;
    state.projectionMatrix(2 * s + 1, o + Acts::eBoundLoc1) = 1.;
  }
  // a dense correlation term, as in a Kalman fit
  Acts::DynamicMatrix b(nPar, nPar);
  for (std::size_t i = 0; i < nPar; ++i) {
    for (std::size_t j = 0; j < nPar; ++j) {
      b(i, j) = std::sin(1. + i + 2. * j);
    }
  }
  state.trackParametersCovariance =
      (state.projectionMatrix.transpose() *
           state.measurementCovariance.inverse() * state.projectionMatrix +
       b.transpose() * b + Acts::DynamicMatrix::Identity(nPar, nPar))
          .inverse();

  state.residual = Acts::DynamicVector::Zero(nMeas);
  state.alignmentToResidualDerivative =
      Acts::DynamicMatrix::Zero(nMeas, state.alignmentDof);
  for (std::size_t i = 0; i < nMeas; ++i) {
    state.residual(i) = 0.01 * (i + 1);
    for (std::size_t a = 0; a < Acts::eAlignmentSize; ++a) {
      state.alignmentToResidualDerivative(i, a) = 0.1 * (a + 1) + 0.01 * i;
    }
  }
  // exactly representable in single precision
  state.residual(0) = static_cast<double>(track + 1);
  return state;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ActsToMilleConcurrentWritingTests)

BOOST_AUTO_TEST_CASE(EveryTrackInItsOwnRecord) {
  const std::size_t nThreads =
      std::max(4u, std::min(16u, std::thread::hardware_concurrency()));
  constexpr std::size_t nTracksPerThread = 2000;
  const std::size_t nTracks = nThreads * nTracksPerThread;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::RectangleBounds>(10., 10.));
  std::vector<ActsAlignment::detail::TrackAlignmentState> states;
  states.reserve(nTracks);
  for (std::size_t track = 0; track < nTracks; ++track) {
    states.push_back(makeState(*surface, track));
  }

  const std::string fname = "ActsToMilleConcurrentWriting.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < nThreads; ++t) {
      threads.emplace_back([&, t]() {
        auto logger = Acts::getDefaultLogger("ActsToMilleConcurrentWriting",
                                             Acts::Logging::WARNING);
        for (std::size_t i = 0; i < nTracksPerThread; ++i) {
          ActsPlugins::ActsToMille::dumpToMille(
              states[t * nTracksPerThread + i], *out, false, *logger);
        }
      });
    }
    for (auto& thread : threads) {
      thread.join();
    }
  }

  // every record holds exactly one track: its surface measurements and one
  // pseudo-measurement per track parameter, the track number first
  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));
  std::size_t nRecords = 0;
  std::size_t nBadRecords = 0;
  std::set<std::size_t> tracksSeen;
  while (true) {
    Mille::MilleDecoder decoder;
    std::vector<Mille::MilleMeasurement> entries;
    const auto res = decoder.decode(*reader, entries);
    if (res == Mille::MilleDecoder::ReadResult::atEof) {
      break;
    }
    BOOST_REQUIRE(res == Mille::MilleDecoder::ReadResult::OK);
    ++nRecords;
    if (entries.size() != nMeas + nPar) {
      ++nBadRecords;
      continue;
    }
    tracksSeen.insert(static_cast<std::size_t>(entries.front().measurement));
  }
  BOOST_CHECK_EQUAL(nBadRecords, 0u);
  BOOST_CHECK_EQUAL(nRecords, nTracks);
  BOOST_CHECK_EQUAL(tracksSeen.size(), nTracks);
}

BOOST_AUTO_TEST_SUITE_END()
