// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Local Trajectory Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <iostream>

#include "Acts/EventData/LocalTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"

using std::cout;
using std::endl;

namespace Acts {
namespace Test {

  using Parameters = TrackParametersBase::ParVector_t;
  using Covariance = TrackParametersBase::CovMatrix_t;

  CurvilinearParameters
  make_params()
  {
    // generate arbitrary positive, definite matrix
    Covariance rnd = Covariance::Random();
    Covariance cov = rnd.transpose() * rnd;
    return {std::make_unique<Covariance>(cov),
            Vector3D(0, 0, 1),
            Vector3D(100, 1000, 400),
            -1};
  }

  BOOST_AUTO_TEST_CASE(local_trajectory_build)
  {
    LocalTrajectory lt;

    // construct trajectory w/ multiple components
    auto i0 = lt.addPoint(make_params());
    // trajectory bifurcates here into multiple hypotheses
    auto i1a = lt.addPoint(make_params(), i0);
    auto i1b = lt.addPoint(make_params(), i0);
    auto i2a = lt.addPoint(make_params(), i1a);
    auto i2b = lt.addPoint(make_params(), i1b);

    // print each trajectory component
    auto print = [](const LocalTrajectoryPoint& p) {
      cout << "  point " << p.index() << endl;
    };
    cout << "trajectory starting at " << i2a << endl;
    lt.traverseBackward(i2a, print);
    cout << "trajectory starting at " << i2b << endl;
    lt.traverseBackward(i2b, print);

    // print full/effective parameters/covariance for an example point
    const auto& p = lt.getPoint(i1b);
    cout << "data for point " << p.index() << endl;
    cout << p.fullParameters().transpose() << endl;
    cout << p.fullCovariance() << endl;
    cout << p.parameters().transpose() << endl;
    cout << p.covariance() << endl;
    cout << "has measurement " << p.hasMeasurement() << endl;
  }

}  // namespace Test
}  // namespace Acts
