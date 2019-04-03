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

    auto i0 = lt.addPoint(make_params());
    auto i1 = lt.addPoint(make_params());
    auto p0 = lt.getPoint(i0);
    auto p1 = lt.getPoint(i1);

    for (const auto& p : {p0, p1}) {
      cout << "point" << endl;
      cout << p.parameters().transpose() << endl;
      cout << p.covariance() << endl;
      cout << p.fullParameters().transpose() << endl;
      cout << p.fullCovariance() << endl;
      cout << "has measurement " << p1.hasMeasurement() << endl;
    }
  }

}  // namespace Test
}  // namespace Acts
