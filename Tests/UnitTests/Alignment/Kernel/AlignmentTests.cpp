// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsTests/CommonHelpers/AlignmentHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <string>

using namespace Acts;
using namespace ActsTests;
using namespace ActsTests::AlignmentUtils;
using namespace Acts::detail::Test;
using namespace Acts::UnitConstants;
///
/// @brief Unit test for KF-based alignment algorithm
///
BOOST_AUTO_TEST_CASE(ZeroFieldKalmanAlignment) {
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

  const auto& alignState = evaluateRes.value();
  CHECK_CLOSE_ABS(alignState.chi2 / alignState.alignmentDof, 0.5, 1);

  // Check the dimensions
  BOOST_CHECK_EQUAL(alignState.measurementDim, 12);
  BOOST_CHECK_EQUAL(alignState.trackParametersDim, 36);
  // Check the alignment dof
  BOOST_CHECK_EQUAL(alignState.alignmentDof, 30);
  BOOST_CHECK_EQUAL(alignState.alignedSurfaces.size(), 5);
  // Check the measurements covariance
  BOOST_CHECK_EQUAL(alignState.measurementCovariance.rows(), 12);
  const SquareMatrix2 measCov =
      alignState.measurementCovariance.block<2, 2>(2, 2);
  SquareMatrix2 cov2D;
  cov2D << 30_um * 30_um, 0, 0, 50_um * 50_um;
  CHECK_CLOSE_ABS(measCov, cov2D, 1e-10);
  // Check the track parameters covariance matrix. Its rows/columns scales
  // with the number of measurement states
  BOOST_CHECK_EQUAL(alignState.trackParametersCovariance.rows(), 36);
  // Check the projection matrix
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.cols(), 36);
  const Matrix<2, 6> proj = alignState.projectionMatrix.block<2, 6>(0, 0);
  const Matrix<2, 6> refProj = Matrix<2, 6>::Identity();
  CHECK_CLOSE_ABS(proj, refProj, 1e-10);
  // Check the residual
  BOOST_CHECK_EQUAL(alignState.residual.size(), 12);
  // Check the residual covariance
  BOOST_CHECK_EQUAL(alignState.residualCovariance.rows(), 12);
  // Check the alignment to residual derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.cols(), 30);
  // Check the chi2 derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2Derivative.size(), 30);
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2SecondDerivative.rows(), 30);

  // Test the align method
  std::vector<std::vector<TestSourceLink>> trajCollection;
  trajCollection.reserve(10);
  std::vector<BoundTrackParameters> sParametersCollection;
  sParametersCollection.reserve(10);
  for (const auto& traj : trajectories) {
    trajCollection.push_back(traj.sourceLinks);
    sParametersCollection.push_back(*traj.startParameters);
  }
  auto alignRes =
      alignZero.align(trajCollection, sParametersCollection, alignOptions);

  // BOOST_CHECK(alignRes.ok());
}
