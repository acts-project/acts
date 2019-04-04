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

#include "Acts/EventData/MultiTrajectory.hpp"
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

  BOOST_AUTO_TEST_CASE(multitrajectory_build)
  {
    MultiTrajectory t;

    // construct trajectory w/ multiple components
    auto i0 = t.addPoint(make_params());
    // trajectory bifurcates here into multiple hypotheses
    auto i1a = t.addPoint(make_params(), i0);
    auto i1b = t.addPoint(make_params(), i0);
    auto i2a = t.addPoint(make_params(), i1a);
    auto i2b = t.addPoint(make_params(), i1b);

    // print each trajectory component
    auto print = [](auto p) {
      cout << "  point " << p.index() << endl;
      cout << "     params " << p.parameters().transpose() << endl;
    };
    cout << "trajectory starting at " << i2a << endl;
    t.visitBackwards(i2a, print);
    cout << "trajectory starting at " << i2b << endl;
    t.visitBackwards(i2b, print);

    // modify elements of the trajectory
    t.applyBackwards(i2b, [](auto p) { p.parameters().setRandom(); });
    cout << "modified trajectory starting at " << i2b << endl;
    t.visitBackwards(i2b, print);

    // print full/effective parameters/covariance for an example point
    const auto& p = t.getPoint(i1b);
    cout << "data for point " << p.index() << endl;
    cout << p.parameters().transpose() << endl;
    cout << p.covariance() << endl;
    cout << "has measurement " << p.hasMeasurement() << endl;
  }

}  // namespace Test
}  // namespace Acts
