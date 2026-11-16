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
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsTests/CommonHelpers/AlignmentHelpers.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <string>
#include <utility>
#include <vector>

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

  // Record proposed updates without moving geometry: each comparison starts
  // from the same detector and fitted tracks.
  using ProposedUpdate = std::pair<SurfacePlacementBase*, Transform3>;
  std::vector<ProposedUpdate> proposedUpdates;
  ActsAlignment::AlignedTransformUpdater voidAlignUpdater =
      [&](SurfacePlacementBase* element, const GeometryContext& /*gctx*/,
          const Transform3& transform) {
        proposedUpdates.emplace_back(element, transform);
        return true;
      };

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
  const auto checkCompatibility =
      [&](const ActsAlignment::AlignmentResult& actual,
          const ActsAlignment::AlignmentResult& expected,
          const std::vector<ProposedUpdate>& expectedUpdates) {
        BOOST_CHECK(actual.idxedAlignSurfaces == expected.idxedAlignSurfaces);
        BOOST_CHECK_EQUAL(actual.alignmentDof, expected.alignmentDof);
        BOOST_CHECK_EQUAL(actual.measurementDim, expected.measurementDim);
        BOOST_CHECK_EQUAL(actual.numTracks, expected.numTracks);
        BOOST_REQUIRE_EQUAL(actual.result.ok(), expected.result.ok());
        if (!expected.result.ok()) {
          BOOST_CHECK_EQUAL(actual.result.error(), expected.result.error());
        }
        CHECK_CLOSE_ABS(actual.chi2, expected.chi2, 1e-12);
        CHECK_CLOSE_ABS(actual.averageChi2ONdf, expected.averageChi2ONdf,
                        1e-12);
        CHECK_CLOSE_ABS(actual.deltaChi2, expected.deltaChi2, 1e-12);
        CHECK_CLOSE_ABS(actual.sumChi2Derivative, expected.sumChi2Derivative,
                        1e-12);
        CHECK_CLOSE_ABS(actual.sumChi2SecondDerivative,
                        expected.sumChi2SecondDerivative, 1e-12);
        CHECK_CLOSE_ABS(actual.deltaAlignmentParameters,
                        expected.deltaAlignmentParameters, 1e-12);
        CHECK_CLOSE_ABS(actual.alignmentCovariance,
                        expected.alignmentCovariance, 1e-12);
        BOOST_REQUIRE_EQUAL(actual.alignedParameters.size(),
                            expected.alignedParameters.size());
        for (const auto& [element, transform] : expected.alignedParameters) {
          const auto it = actual.alignedParameters.find(element);
          BOOST_REQUIRE(it != actual.alignedParameters.end());
          CHECK_CLOSE_ABS(it->second.matrix(), transform.matrix(), 1e-12);
        }
        BOOST_REQUIRE_EQUAL(proposedUpdates.size(), expectedUpdates.size());
        for (std::size_t i = 0; i < expectedUpdates.size(); ++i) {
          BOOST_CHECK(proposedUpdates[i].first == expectedUpdates[i].first);
          CHECK_CLOSE_ABS(proposedUpdates[i].second.matrix(),
                          expectedUpdates[i].second.matrix(), 1e-12);
        }
      };

  // Cover both immediate convergence and one update followed by convergence.
  // The recording updater leaves chi2 unchanged, so the latter converges when
  // the second iteration compares it with the first.
  for (bool exerciseUpdate : {false, true}) {
    BOOST_TEST_CONTEXT("exerciseUpdate=" << exerciseUpdate) {
      alignOptions.maxIterations = exerciseUpdate ? 2 : 1;
      alignOptions.averageChi2ONdfCutOff = exerciseUpdate ? 0. : 0.5;
      alignOptions.deltaAverageChi2ONdfCutOff = {1, 0.};
      alignOptions.alignedStructures.clear();
      alignOptions.activeStructures.clear();
      proposedUpdates.clear();
      auto alignRes =
          alignZero.align(trajCollection, sParametersCollection, alignOptions);
      BOOST_REQUIRE(alignRes.ok());
      BOOST_REQUIRE(alignRes->result.ok());
      BOOST_REQUIRE_GT(alignRes->measurementDim, 0u);
      BOOST_REQUIRE_EQUAL(alignRes->deltaAlignmentParameters.size(),
                          alignRes->alignmentDof);
      BOOST_REQUIRE(alignRes->deltaAlignmentParameters.allFinite());
      BOOST_REQUIRE(alignRes->alignmentCovariance.allFinite());
      const auto expectedUpdates = proposedUpdates;
      const auto expectedUpdateCount =
          exerciseUpdate ? alignOptions.alignedDetElements.size() : 0u;
      BOOST_REQUIRE_EQUAL(expectedUpdates.size(), expectedUpdateCount);
      BOOST_REQUIRE_EQUAL(alignRes->alignedParameters.size(),
                          expectedUpdateCount);

      // Registration alone does not activate group solving, even with All and
      // no structure frame provider.
      auto metadata = std::make_shared<ActsAlignment::AlignableStructure>(
          GeometryIdentifier());
      metadata->alignmentMask() = ActsAlignment::AlignmentMask::All;
      metadata->addSurface(
          detector.detectorStore.front()->surface().getSharedPtr());
      alignOptions.alignedStructures = {metadata};
      proposedUpdates.clear();
      auto metadataResult =
          alignZero.align(trajCollection, sParametersCollection, alignOptions);
      BOOST_REQUIRE(metadataResult.ok());
      checkCompatibility(*metadataResult, *alignRes, expectedUpdates);

      // Explicit None preserves the same module results and update proposals.
      metadata->alignmentMask() = ActsAlignment::AlignmentMask::None;
      alignOptions.activeStructures = {metadata};
      proposedUpdates.clear();
      auto inactiveResult =
          alignZero.align(trajCollection, sParametersCollection, alignOptions);
      BOOST_REQUIRE(inactiveResult.ok());
      checkCompatibility(*inactiveResult, *alignRes, expectedUpdates);
    }
  }
}

namespace {
struct SelectionTestOptions {
  GeometryContext geoContext = GeometryContext::dangerouslyDefaultConstruct();
  const Surface* referenceSurface = nullptr;
};

struct SelectionTestTrack {
  std::size_t tipIndex() const { return 0; }
};

struct CountingSelectionFitter {
  std::size_t* calls;

  template <typename Iterator, typename Parameters, typename Options,
            typename Tracks>
  Result<SelectionTestTrack> fit(Iterator /*begin*/, Iterator /*end*/,
                                 const Parameters& /*parameters*/,
                                 const Options& /*options*/,
                                 Tracks& /*tracks*/) const {
    ++*calls;
    return std::make_error_code(std::errc::operation_canceled);
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(StructureSolveGuardPrecedesFittingAndUpdates) {
  using namespace ActsAlignment;
  auto element = std::make_shared<DetectorElementStub>(
      Transform3::Identity(),
      std::make_shared<const RectangleBounds>(10_mm, 10_mm), 1_um);
  auto structure = std::make_shared<AlignableStructure>(GeometryIdentifier(),
                                                        Transform3::Identity());
  structure->addSurface(element->surface().getSharedPtr());
  structure->alignmentMask() = AlignmentMask::All;
  std::size_t fitCalls = 0;
  std::size_t updateCalls = 0;
  Alignment engine(CountingSelectionFitter{&fitCalls},
                   getDummyLogger().clone());
  AlignedTransformUpdater updater = [&](auto*, const auto&, const auto&) {
    ++updateCalls;
    return true;
  };
  AlignmentOptions options(SelectionTestOptions{}, updater, {}, 0.5, {5, 0.01},
                           1, {}, {structure}, {structure});
  const std::vector<std::vector<SourceLink>> tracks(1);
  const std::vector<BoundTrackParameters> parameters{
      BoundTrackParameters::createCurvilinear(Vector4::Zero(), Vector3::UnitZ(),
                                              1., std::nullopt,
                                              ParticleHypothesis::pion())};
  auto result = engine.align(tracks, parameters, options);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    AlignmentError::UnsupportedStructureAlignment);
  BOOST_CHECK_EQUAL(fitCalls, 0u);
  BOOST_CHECK_EQUAL(updateCalls, 0u);

  options.activeStructures = {nullptr};
  auto invalid = engine.align(tracks, parameters, options);
  BOOST_REQUIRE(!invalid.ok());
  BOOST_CHECK_EQUAL(invalid.error(), AlignmentError::InvalidAlignmentSelection);

  AlignableStructure cycle{GeometryIdentifier()};
  auto alias = std::shared_ptr<AlignableStructure>(&cycle, [](auto*) {});
  cycle.addChild(alias);
  options.alignedStructures = {alias};
  options.activeStructures.clear();
  auto cyclic = engine.align(tracks, parameters, options);
  BOOST_REQUIRE(!cyclic.ok());
  BOOST_CHECK_EQUAL(cyclic.error(), AlignmentError::HierarchyValidationFailure);
  BOOST_CHECK_EQUAL(fitCalls, 0u);
  BOOST_CHECK_EQUAL(updateCalls, 0u);
}
