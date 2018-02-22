// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * @file GainMatrixTests.cpp
 */

// STL include(s)
#include <memory>

// Boost include(s)
#define BOOST_TEST_MODULE Measurement Tests
#include <boost/test/included/unit_test.hpp>

// ATS include(s)
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Fitter/KalmanUpdator.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {
  template <ParID_t... params>
  using Measurement_t = Measurement<unsigned long int, params...>;

  BOOST_AUTO_TEST_CASE(gain_matrix_updator)
  {
    // make dummy measurement
    CylinderSurface   cylinder(nullptr, 3, 10);
    ActsSymMatrixD<2> cov;
    cov << 0.04, 0, 0, 0.1;
    FittableMeasurement<unsigned long int> m
        = Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1>(
            cylinder, 0, std::move(cov), -0.1, 0.45);

    // make dummy track parameter
    ActsSymMatrixD<Acts::NGlobalPars> covTrk;
    covTrk << 0.08, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
    ActsVectorD<Acts::NGlobalPars> parValues;
    parValues << 0.3, 0.5, 0.5 * M_PI, 0.3 * M_PI, 0.01;
    BoundParameters pars(
        std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
        parValues,
        cylinder);

    GainMatrixUpdator                      g;
    std::unique_ptr<const BoundParameters> filtered(g(m, pars));

    std::cout << pars << std::endl;
    if (filtered) std::cout << *filtered << std::endl;
  }
}  // end of namespace Test
}  // end of namespace Acts
