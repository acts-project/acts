// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE TrackState Tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

  using Identifier = unsigned long int;

  template <ParID_t... params>
  using MeasurementType = Measurement<Identifier, params...>;

  template <ParID_t... params>
  using MeasuredStateType
      = MeasuredTrackState<Identifier, TrackParameters, params...>;
  using ParametricStateType = ParametricTrackState<Identifier, TrackParameters>;
  using VariantStateType    = VariantTrackState<Identifier, TrackParameters>;
  ///
  /// @brief Unit test for creation of Measurement object
  ///
  BOOST_AUTO_TEST_CASE(track_state_initialization)
  {
    // The plane surface
    PlaneSurface plane(Vector3D{0., 0., 0.}, Vector3D{1., 0., 0.});

    // Construct the 1D measurement
    ActsSymMatrixD<1> cov1D;
    cov1D << 0.04;
    MeasurementType<ParDef::eLOC_0> m1D(plane, 0, std::move(cov1D), 0.02);
    // Construct the 2D measurement
    ActsSymMatrixD<2> cov2D;
    cov2D << 0.04, 0., 0.09, 0.;
    MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m2D(
        plane, 0, std::move(cov2D), 0.02, 0.03);

    // The 1D track state from the measurement
    VariantStateType mts1D = MeasuredStateType<ParDef::eLOC_0>(m1D);
    // The 2D track state from the measurement
    VariantStateType mts2D
        = MeasuredStateType<ParDef::eLOC_0, ParDef::eLOC_1>(m2D);

    // Construct the parameter
    std::array<double, 5> pars_array = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    // constructor from parameter vector
    auto ataPlane
        = std::make_unique<const BoundParameters>(nullptr, pars, plane);
    // The parametric track state from the parameters
    VariantStateType pts = ParametricStateType(std::move(ataPlane));

    // The parametric track state from the surface only
    VariantStateType sts = ParametricStateType(plane);
  }

}  // namespace Test
}  // namespace Acts
