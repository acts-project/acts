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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "Mille/MilleFactory.h"

// Reading several records into the same TrackAlignmentState, as
// ActsSolverFromMille does. Tracks do not all cross the same surfaces
// (holes, tracks leaving the detector): every record must be unpacked on its
// own, with the derivative matrix sized and indexed for its own surfaces.

namespace {

constexpr std::size_t nSurfaces = 6;

/// one measurement: residual, the surface it is on, the alignment labels
/// (dof 0..5) with non-zero derivative
struct Hit {
  float residual;
  std::size_t surface;
  std::vector<std::size_t> dofs = {0, 1, 2, 3, 4, 5};
};

int label(std::size_t surface, std::size_t dof) {
  return static_cast<int>(surface * Acts::eAlignmentSize + dof + 1);
}
double derivative(std::size_t surface, std::size_t dof) {
  return 0.25 + 0.5 * surface + 0.125 * dof;
}

void writeRecord(Mille::MilleRecord& record, const std::vector<Hit>& hits) {
  unsigned int local = 1;
  for (const Hit& hit : hits) {
    std::vector<int> labels;
    std::vector<double> derivatives;
    for (std::size_t dof : hit.dofs) {
      labels.push_back(label(hit.surface, dof));
      derivatives.push_back(derivative(hit.surface, dof));
    }
    // one local parameter per measurement, so that the local system is
    // well defined
    record.addData(hit.residual, 0.1f, std::vector<unsigned int>{local},
                   std::vector<double>{1.0}, labels, derivatives);
    ++local;
  }
  record.writeRecord();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ActsMilleReadRecordSequenceTests)

BOOST_AUTO_TEST_CASE(ReuseStateAcrossRecords) {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("ReuseStateAcrossRecords", Acts::Logging::INFO));
  // the alignable surfaces and their indices in the geometry
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  std::unordered_map<const Acts::Surface*, std::size_t> idxedAlignSurfaces;
  for (std::size_t s = 0; s < nSurfaces; ++s) {
    surfaces.push_back(Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::RectangleBounds>(10., 10.)));
    idxedAlignSurfaces.emplace(surfaces.back().get(), s);
  }

  // record sequence: a full track first, then tracks with fewer surfaces
  const std::vector<std::vector<Hit>> records = {
      // full track
      {{0.01f, 0}, {0.02f, 1}, {0.03f, 2}, {0.04f, 3}, {0.05f, 4}, {0.06f, 5}},
      // track on the last surface only
      {{0.07f, 5}},
      // track with a hole on surface 2
      {{0.01f, 0}, {0.02f, 1}, {0.04f, 3}, {0.05f, 4}, {0.06f, 5}},
      // track leaving after surface 1
      {{0.01f, 0}, {0.02f, 1}},
      // surface 4 seen only through dof 2..5: its first label is missing
      {{0.03f, 3}, {0.05f, 4, {2, 3, 4, 5}}},
      // a hit with a residual of exactly zero is still a hit
      {{0.f, 1}, {0.02f, 2}},
  };

  const std::string fname = "ActsMilleReadRecordSequence.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    for (const auto& hits : records) {
      writeRecord(*out, hits);
    }
  }

  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));

  // one state for all records, as in ActsSolverFromMille
  ActsAlignment::detail::TrackAlignmentState state;
  for (const auto& hits : records) {
    BOOST_REQUIRE(ActsPlugins::ActsToMille::unpackMilleRecord(
                      *reader, state, idxedAlignSurfaces, logger()) ==
                  Mille::MilleDecoder::ReadResult::OK);

    std::set<std::size_t> onTrack;
    for (const Hit& hit : hits) {
      onTrack.insert(hit.surface);
    }
    // exactly the surfaces of this record, numbered 0..k-1 in label order
    BOOST_CHECK_EQUAL(state.alignedSurfaces.size(), onTrack.size());
    BOOST_CHECK_EQUAL(state.alignmentDof,
                      Acts::eAlignmentSize * onTrack.size());
    BOOST_CHECK_EQUAL(state.measurementDim, hits.size());
    BOOST_CHECK_EQUAL(state.alignmentToResidualDerivative.rows(),
                      static_cast<long>(hits.size()));
    BOOST_CHECK_EQUAL(state.alignmentToResidualDerivative.cols(),
                      static_cast<long>(state.alignmentDof));
    std::map<std::size_t, std::size_t> internalOf;
    std::size_t expectedInternal = 0;
    for (std::size_t s : onTrack) {
      internalOf[s] = expectedInternal++;
    }
    for (const auto& [surface, indices] : state.alignedSurfaces) {
      const auto [globalIndex, internalIndex] = indices;
      BOOST_CHECK_EQUAL(globalIndex, idxedAlignSurfaces.at(surface));
      BOOST_REQUIRE(internalOf.contains(globalIndex));
      BOOST_CHECK_EQUAL(internalIndex, internalOf.at(globalIndex));
    }
    // every derivative in the column of its label, zero elsewhere
    for (std::size_t iMeas = 0; iMeas < hits.size(); ++iMeas) {
      const Hit& hit = hits[iMeas];
      for (std::size_t s : onTrack) {
        for (std::size_t dof = 0; dof < Acts::eAlignmentSize; ++dof) {
          const bool written =
              s == hit.surface &&
              std::ranges::find(hit.dofs, dof) != hit.dofs.end();
          const double value = state.alignmentToResidualDerivative(
              iMeas, internalOf.at(s) * Acts::eAlignmentSize + dof);
          if (written) {
            BOOST_CHECK_CLOSE(value, derivative(s, dof), 1e-4);
          } else {
            BOOST_CHECK_EQUAL(value, 0.);
          }
        }
      }
    }
  }
  BOOST_CHECK(ActsPlugins::ActsToMille::unpackMilleRecord(
                  *reader, state, idxedAlignSurfaces, logger()) ==
              Mille::MilleDecoder::ReadResult::atEof);
}

BOOST_AUTO_TEST_CASE(UnknownSurfaceIsReadError) {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("UnknownSurfaceIsReadError", Acts::Logging::INFO));
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  std::unordered_map<const Acts::Surface*, std::size_t> idxedAlignSurfaces;
  for (std::size_t s = 0; s < 2; ++s) {
    surfaces.push_back(Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::RectangleBounds>(10., 10.)));
    idxedAlignSurfaces.emplace(surfaces.back().get(), s);
  }

  const std::string fname = "ActsMilleReadUnknownSurface.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    writeRecord(*out, {{0.01f, 0}, {0.02f, 1}});
    // surface index 4 is not in the indexed list
    writeRecord(*out, {{0.01f, 0}, {0.02f, 4}});
  }
  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));

#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  // the read error below is expected: do not fail on it
  Acts::Logging::ScopedFailureThreshold threshold{Acts::Logging::FATAL};
#endif
  ActsAlignment::detail::TrackAlignmentState state;
  BOOST_REQUIRE(ActsPlugins::ActsToMille::unpackMilleRecord(
                    *reader, state, idxedAlignSurfaces, logger()) ==
                Mille::MilleDecoder::ReadResult::OK);
  BOOST_CHECK_EQUAL(state.alignmentDof, 2 * Acts::eAlignmentSize);
  BOOST_CHECK(ActsPlugins::ActsToMille::unpackMilleRecord(
                  *reader, state, idxedAlignSurfaces, logger()) ==
              Mille::MilleDecoder::ReadResult::error);
  // a failed read leaves the state untouched
  BOOST_CHECK_EQUAL(state.alignmentDof, 2 * Acts::eAlignmentSize);
  BOOST_CHECK_EQUAL(state.alignedSurfaces.size(), 2u);
}

BOOST_AUTO_TEST_CASE(GlobalsWithoutLocalsIsReadError) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("GlobalsWithoutLocalsIsReadError",
                                           Acts::Logging::INFO));
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  std::unordered_map<const Acts::Surface*, std::size_t> idxedAlignSurfaces;
  for (std::size_t s = 0; s < 2; ++s) {
    surfaces.push_back(Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::RectangleBounds>(10., 10.)));
    idxedAlignSurfaces.emplace(surfaces.back().get(), s);
  }

  const std::string fname = "ActsMilleReadGlobalsWithoutLocals.dat";
  {
    std::unique_ptr<Mille::MilleRecord> out = Mille::spawnMilleRecord(fname);
    BOOST_REQUIRE(out != nullptr);
    writeRecord(*out, {{0.01f, 0}, {0.02f, 1}});
    // second record: the hit on surface 1 has alignment derivatives but its
    // only local derivative is zero, which Mille does not write
    out->addData(0.01f, 0.1f, std::vector<unsigned int>{1u},
                 std::vector<double>{1.0}, std::vector<int>{label(0, 0)},
                 std::vector<double>{derivative(0, 0)});
    out->addData(0.02f, 0.1f, std::vector<unsigned int>{2u},
                 std::vector<double>{0.0}, std::vector<int>{label(1, 0)},
                 std::vector<double>{derivative(1, 0)});
    out->writeRecord();
  }
  auto reader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(reader != nullptr);
  BOOST_REQUIRE(reader->open(fname));

#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  // the read error below is expected: do not fail on it
  Acts::Logging::ScopedFailureThreshold threshold{Acts::Logging::FATAL};
#endif
  ActsAlignment::detail::TrackAlignmentState state;
  BOOST_REQUIRE(ActsPlugins::ActsToMille::unpackMilleRecord(
                    *reader, state, idxedAlignSurfaces, logger()) ==
                Mille::MilleDecoder::ReadResult::OK);
  BOOST_CHECK(ActsPlugins::ActsToMille::unpackMilleRecord(
                  *reader, state, idxedAlignSurfaces, logger()) ==
              Mille::MilleDecoder::ReadResult::error);
  // a failed read leaves the state untouched
  BOOST_CHECK_EQUAL(state.alignmentDof, 2 * Acts::eAlignmentSize);
}

BOOST_AUTO_TEST_SUITE_END()
