// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE ZScanVertexFinder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Vertexing/MultiAdaptiveVertexFitter.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

/// @brief Unit test for MultiAdaptiveVertexFitter
///
BOOST_AUTO_TEST_CASE(multi_adaptive_vertex_fitter_test) {
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  Propagator<EigenStepper<ConstantBField>> propagator(stepper);

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>::Config
      config(bField, propagator);

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>::State
      state;

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>
      fitter(config);

  Vertex<BoundParameters> vtx;

  auto res1 = fitter.addVertexToFit(state, vtx, fitterOptions);

  auto res2 = fitter.fit(state, fitterOptions);
}

}  // namespace Test
}  // namespace Acts