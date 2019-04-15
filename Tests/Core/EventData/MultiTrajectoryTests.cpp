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
#include "Acts/EventData/TrackState.hpp"

using std::cout;
using std::endl;

namespace Acts {
namespace Test {

  using SourceLink = unsigned long long;
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

  TrackState<SourceLink, CurvilinearParameters>
  make_trackstate()
  {
    return {make_params()};
  }

  BOOST_AUTO_TEST_CASE(multitrajectory_build)
  {
    MultiTrajectory<SourceLink> t;

    // construct trajectory w/ multiple components
    auto i0 = t.addTrackState(make_trackstate());
    // trajectory bifurcates here into multiple hypotheses
    auto i1a = t.addTrackState(make_trackstate(), i0);
    auto i1b = t.addTrackState(make_trackstate(), i0);
    auto i2a = t.addTrackState(make_trackstate(), i1a);
    auto i2b = t.addTrackState(make_trackstate(), i1b);

    // print each trajectory component
    auto print = [](auto p) {
      cout << "  point " << p.index() << endl;
      cout << "     params " << p.predicted().transpose() << endl;
    };
    cout << "trajectory starting at " << i2a << endl;
    t.visitBackwards(i2a, print);
    cout << "trajectory starting at " << i2b << endl;
    t.visitBackwards(i2b, print);

    // modify elements of the trajectory
    t.applyBackwards(i2b, [](auto p) { p.predicted().setRandom(); });
    cout << "modified trajectory starting at " << i2b << endl;
    t.visitBackwards(i2b, print);

    // print full/effective parameters/covariance for an example point
    const auto& p = t.getTrackState(i1b);
    cout << "data for point " << p.index() << endl;
    cout << p.predicted().transpose() << endl;
    cout << p.predictedCovariance() << endl;
    cout << "has uncalibrated " << p.hasUncalibrated() << endl;
  }

}  // namespace Test
}  // namespace Acts
