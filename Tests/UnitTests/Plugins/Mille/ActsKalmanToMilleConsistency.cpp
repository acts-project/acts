// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"
#include "ActsPlugins/Mille/Helpers.hpp"
#include "ActsTests/CommonHelpers/AlignmentHelpers.hpp"

#include <string>

#include <Mille/MilleDecoder.h>
#include <Mille/MilleFactory.h>

using namespace Acts;
using namespace ActsTests;
using namespace ActsTests::AlignmentUtils;
using namespace Acts::detail::Test;
using namespace Acts::UnitConstants;

BOOST_AUTO_TEST_CASE(ZeroFieldKalmanToMille) {
  /// Part 1): Build some test data.
  /// This is shamelessly "borrowed" from
  /// the existing alignment unit test.

  aliTestUtils utils;
  // Build detector
  TelescopeDetector detector(utils.geoCtx);
  const auto geometry = detector();

  // reconstruction propagator and fitter
  auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
  const auto kfZeroPropagator =
      makeConstantFieldPropagator(geometry, 0_T, std::move(kfLogger));
  auto kfZero = KalmanFitterType(kfZeroPropagator);

  // alignment
  auto alignLogger = getDefaultLogger("Alignment", Logging::INFO);
  const auto alignZero =
      ActsAlignment::Alignment(std::move(kfZero), std::move(alignLogger));

  // Create 10 trajectories
  const auto& trajectories = createTrajectories(geometry, 10, utils);

  // Construct the KalmanFitter options

  auto extensions = getExtensions(utils);
  TestSourceLink::SurfaceAccessor surfaceAccessor{*geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
  KalmanFitterOptions kfOptions(
      utils.geoCtx, utils.magCtx, utils.calCtx, extensions,
      PropagatorPlainOptions(utils.geoCtx, utils.magCtx));

  // Construct a non-updating alignment updater
  ActsAlignment::AlignedTransformUpdater voidAlignUpdater =
      [](SurfacePlacementBase* /*element*/, const GeometryContext& /*gctx*/,
         const Transform3& /*transform*/) { return true; };

  // Construct the alignment options
  ActsAlignment::AlignmentOptions<KalmanFitterOptions<VectorMultiTrajectory>>
      alignOptions(kfOptions, voidAlignUpdater);
  alignOptions.maxIterations = 1;

  // Set the surfaces to be aligned (fix the layer 8)
  unsigned int iSurface = 0;
  std::unordered_map<const Surface*, std::size_t> idxedAlignSurfaces;
  // Loop over the detector elements
  for (auto& det : detector.detectorStore) {
    const auto& surface = det->surface();
    if (surface.geometryId().layer() != 8) {
      alignOptions.alignedDetElements.push_back(det.get());
      idxedAlignSurfaces.emplace(&surface, iSurface);
      iSurface++;
    }
  }

  // Test the method to evaluate alignment state for a single track
  const auto& inputTraj = trajectories.front();
  kfOptions.referenceSurface = &(*inputTraj.startParameters).referenceSurface();

  auto evaluateRes = alignZero.evaluateTrackAlignmentState(
      kfOptions.geoContext, inputTraj.sourceLinks, *inputTraj.startParameters,
      kfOptions, idxedAlignSurfaces, ActsAlignment::AlignmentMask::All);
  BOOST_CHECK(evaluateRes.ok());

  /// Part 2: Now dump the alignment state to Mille.

  auto milleRecord = Mille::spawnMilleRecord("myRecord.root", true);

  const auto& alignState = evaluateRes.value();
  ActsPlugins::ActsToMille::dumpToMille(alignState, *milleRecord);

  // trigger file close by destroying the Mille record
  milleRecord->flushOutputFile();
  milleRecord.reset();

  /// Part 3: We read back the record from the binary file,
  /// convert it to an alignment state, and compare to
  /// what we started from

  // read back the record
  auto milleReader = Mille::spawnMilleReader("myRecord.root");
  BOOST_CHECK(milleReader->open("myRecord.root"));

  // Decode it back into a TrackAlignmentState

  ActsAlignment::detail::TrackAlignmentState millePedeState;
  // we need to externally supply the alignment parameter indexing logic
  millePedeState.alignedSurfaces = alignState.alignedSurfaces;
  BOOST_CHECK(ActsPlugins::ActsToMille::unpackMilleRecord(
                  *milleReader, millePedeState, idxedAlignSurfaces) ==
              Mille::MilleDecoder::ReadResult::OK);

  // now compare the results!

  BOOST_CHECK_EQUAL(millePedeState.alignmentDof, alignState.alignmentDof);
  BOOST_CHECK_EQUAL(millePedeState.measurementDim, alignState.measurementDim);
  BOOST_CHECK_EQUAL(millePedeState.trackParametersDim,
                    alignState.trackParametersDim);

  BOOST_CHECK_EQUAL(millePedeState.residual, alignState.residual);
  BOOST_CHECK_EQUAL(millePedeState.chi2, alignState.chi2);
  BOOST_CHECK_EQUAL(millePedeState.measurementCovariance,
                    alignState.measurementCovariance);
  BOOST_CHECK_EQUAL(millePedeState.projectionMatrix,
                    alignState.projectionMatrix);

  BOOST_CHECK_EQUAL(millePedeState.alignmentToResidualDerivative,
                    alignState.alignmentToResidualDerivative);

  // for the track parameter covariance, we only compare the regularised version
  // (differences expected in poorly constrained parameters. Also increase the
  // tolerance a bit.
  BOOST_CHECK(millePedeState.trackParametersCovariance.isApprox(
      ActsPlugins::ActsToMille::regulariseCovariance(
          alignState.trackParametersCovariance),
      1.e-6));

  // for anything derived from the TP covariance, we can also only ask for
  // approximate equivalence as a result.
  BOOST_CHECK(millePedeState.residualCovariance.isApprox(
      alignState.residualCovariance, 1.e-6));
  BOOST_CHECK(millePedeState.alignmentToChi2Derivative.isApprox(
      alignState.alignmentToChi2Derivative, 1.e-6));
  BOOST_CHECK(millePedeState.alignmentToChi2SecondDerivative.isApprox(
      alignState.alignmentToChi2SecondDerivative, 1.e-6));
}
