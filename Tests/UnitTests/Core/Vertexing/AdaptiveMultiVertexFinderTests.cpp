// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPoint3dEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator, pOptions, false);
  IPEstimator ip3dEst(ip3dEstCfg);

  std::vector<double> temperatures(1, 3.);
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundParameters, Linearizer>;

  Fitter::Config fitterCfg(ip3dEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator, pOptions);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  // fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = TrackDensityVertexFinder<Fitter, GaussianTrackDensity>;

  SeedFinder seedFinder;

  using IPEstimater = TrackToVertexIPEstimator<BoundParameters, Propagator>;

  IPEstimater::Config ipEstCfg(propagator, pOptions);

  // Create TrackToVertexIPEstimator
  IPEstimater ipEst(ipEstCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), std::move(seedFinder),
                              std::move(ipEst), std::move(linearizer));

  Finder finder(finderConfig);

  std::vector<BoundParameters> tracks;

  VertexFinderOptions<BoundParameters> finderOptions(tgContext, mfContext);

  std::cout << &finderOptions << std::endl;

  Transform3D transform;
  ActsSymMatrixD<3> rotMat;
  rotMat << -0.780869, 0.624695, 3.79565e-17, 0.455842, 0.569803, 0.683763,
      0.427144, 0.53393, -0.729704;
  transform.rotate(rotMat);
  transform.translation() = Vector3D(1, 2, 3);

  if (debugMode) {
    std::cout << "Debug mode." << std::endl;
  }
  std::cout << "translation: " << transform.translation() << std::endl;
  std::cout << "rot: " << transform.rotation() << std::endl;

  auto bla = finder.find(tracks, finderOptions);
}

}  // namespace Test
}  // namespace Acts
