// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Helpers Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <bitset>

using namespace Acts::VectorHelpers;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Utilities)

  BOOST_AUTO_TEST_CASE(bitset_to_matrix_to_bitset)
  {
    Eigen::Matrix<int, 4, 3> mat;
    mat << 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0;

    std::bitset<4 * 3> act = matrixToBitset(mat);
    std::bitset<4 * 3> exp{"101001011010"};

    BOOST_CHECK_EQUAL(exp, act);

    Eigen::Matrix<int, 4, 3> cnv;
    cnv = bitsetToMatrix<decltype(cnv)>(act);

    BOOST_CHECK_EQUAL(mat, cnv);
  }

struct MyStruct {
  double phi() const { return 42; }
};

BOOST_AUTO_TEST_CASE(phi_helper_test) {
  Vector3D v(0, 1, 0);
  CHECK_CLOSE_ABS(phi(v), M_PI / 2., 1e-9);

  // should work with dynamic types as well
  ActsVectorXd v2{3};
  v2 << 0, 1, 0;
  CHECK_CLOSE_ABS(phi(v2), M_PI / 2., 1e-9);

  MyStruct s;
  BOOST_CHECK_EQUAL(phi(s), 42);
}

BOOST_AUTO_TEST_CASE(perp_helper_test) {
  Vector3D v(1, 2, 3);
  CHECK_CLOSE_ABS(perp(v), std::sqrt(1 + 2 * 2), 1e-9);

  // should work with dynamic types as well
  ActsVectorXd v2{3};
  v2 << 1, 2, 3;
  CHECK_CLOSE_ABS(perp(v2), std::sqrt(1 + 2 * 2), 1e-9);
}

BOOST_AUTO_TEST_CASE(theta_eta_test_helper) {
  Vector3D v(1, 2, 3);
  CHECK_CLOSE_ABS(theta(v), 0.640522, 1e-5);
  CHECK_CLOSE_ABS(eta(v), 1.10359, 1e-5);

  // should work with dynamic types as well
  ActsVectorXd v2{3};
  v2 << 1, 2, 3;
  CHECK_CLOSE_ABS(theta(v2), 0.640522, 1e-5);
  CHECK_CLOSE_ABS(eta(v2), 1.10359, 1e-5);
}

BOOST_AUTO_TEST_CASE(cross_test_helper) {
  {
    Vector3D v(1, 2, 3);
    ActsMatrixD<3, 3> mat;
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    ActsMatrixD<3, 3> act = cross(mat, v);
    ActsMatrixD<3, 3> exp;
    exp << -2, -1, 0, 4, 2, 0, -2, -1, 0;

    CHECK_CLOSE_ABS(act, exp, 1e-9);
  }

  // should work with dynamic types as well
  {
    ActsVectorXd v{3};
    v << 1, 2, 3;
    ActsMatrixXd mat{3, 3};
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    ActsMatrixXd act = cross(mat, v);
    ActsMatrixXd exp{3, 3};
    exp << -2, -1, 0, 4, 2, 0, -2, -1, 0;

    BOOST_CHECK(act.isApprox(exp, 1e-9));
  }
}

BOOST_AUTO_TEST_CASE(toString_test_helper) {
  ActsMatrixD<3, 3> mat;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::string out;
  out = toString(mat);
  BOOST_CHECK(out.size() > 0);

  Translation3D trl{Vector3D{1, 2, 3}};
  out = toString(trl);
  BOOST_CHECK(out.size() > 0);

  Transform3D trf;
  trf = trl;
  out = toString(trf);
  BOOST_CHECK(out.size() > 0);
}

BOOST_AUTO_TEST_CASE(shared_vector_helper_test) {
  {
    std::vector<std::shared_ptr<int>> vec;
    vec = {std::make_shared<int>(5), std::make_shared<int>(9),
           std::make_shared<int>(26), std::make_shared<int>(18473)};

    std::vector<int*> unpacked = unpack_shared_vector(vec);

    std::vector<int*> exp = {
        vec[0].get(),
        vec[1].get(),
        vec[2].get(),
        vec[3].get(),
    };

    BOOST_CHECK_EQUAL_COLLECTIONS(unpacked.begin(), unpacked.end(), exp.begin(),
                                  exp.end());
  }

  // same for const
  {
    std::vector<std::shared_ptr<const int>> vec;
    vec = {std::make_shared<const int>(5), std::make_shared<const int>(9),
           std::make_shared<const int>(26), std::make_shared<const int>(18473)};

    std::vector<const int*> unpacked = unpack_shared_vector(vec);

    std::vector<const int*> exp = {
        vec[0].get(),
        vec[1].get(),
        vec[2].get(),
        vec[3].get(),
    };

    BOOST_CHECK_EQUAL_COLLECTIONS(unpacked.begin(), unpacked.end(), exp.begin(),
                                  exp.end());
  }
}

BOOST_AUTO_TEST_CASE(time_helper_test) {
  {
    SpacePointVector spaceTimeVec = {1., 2., 3., 4.};
    double t = time(spaceTimeVec);
    BOOST_CHECK_EQUAL(spaceTimeVec[3], t);

    double newTime = 17;
    time(spaceTimeVec) = newTime;
    BOOST_CHECK_EQUAL(spaceTimeVec[3], newTime);

    BoundVector boundVec;
    boundVec << 1., 2., 3., 4., 5., 6.;
    t = time(boundVec);
    BOOST_CHECK_EQUAL(boundVec[eT], t);

    newTime = 28;
    time(boundVec) = newTime;
    BOOST_CHECK_EQUAL(boundVec[eT], newTime);

    FreeVector freeVec;
    freeVec << 1., 2., 3., 4., 5., 6., 7., 8.;
    t = time(freeVec);
    BOOST_CHECK_EQUAL(freeVec[7], t);

    newTime = 352;
    time(freeVec) = newTime;
    BOOST_CHECK_EQUAL(freeVec[7], newTime);
  }

  // same for const
  {
    // use temp object to initialize const vectors
    SpacePointVector tempSpaceVec = {1., 2., 3., 4.};
    const SpacePointVector spaceTimeVec(tempSpaceVec);
    const double t1 = time(spaceTimeVec);
    BOOST_CHECK_EQUAL(spaceTimeVec[3], t1);

    BoundVector tempBoundVec;
    tempBoundVec << 1., 2., 3., 4., 5., 6.;
    const BoundVector boundVec(tempBoundVec);
    const double t2 = time(boundVec);
    BOOST_CHECK_EQUAL(boundVec[eT], t2);

    FreeVector tempFreeVec;
    tempFreeVec << 1., 2., 3., 4., 5., 6., 7., 8.;
    const FreeVector freeVec(tempFreeVec);
    const double t3 = time(freeVec);
    BOOST_CHECK_EQUAL(freeVec[7], t3);
  }
}

BOOST_AUTO_TEST_CASE(position_helper_test) {
  {
    SpacePointVector spaceTimeVec = {1., 2., 3., 4.};
    Vector3D pos = position(spaceTimeVec);
    BOOST_CHECK_EQUAL(spaceTimeVec.head<3>(), pos);

    Vector3D newPos = {4., 7., 12.};
    position(spaceTimeVec) = newPos;
    BOOST_CHECK_EQUAL(spaceTimeVec.head<3>(), newPos);

    FreeVector freeVec;
    freeVec << 1., 2., 3., 4., 5., 6., 7., 8.;
    pos = position(freeVec);
    BOOST_CHECK_EQUAL(freeVec.head<3>(), pos);

    newPos = {14., 17., 112.};
    position(freeVec) = newPos;
    BOOST_CHECK_EQUAL(freeVec.head<3>(), newPos);
  }

  // same for const
  {
    // use temp object to initialize const vectors
    SpacePointVector tempSpaceVec = {1., 2., 3., 4.};
    const SpacePointVector spaceTimeVec(tempSpaceVec);
    const Vector3D pos1 = position(spaceTimeVec);
    BOOST_CHECK_EQUAL(spaceTimeVec.head<3>(), pos1);

    FreeVector tempFreeVec;
    tempFreeVec << 1., 2., 3., 4., 5., 6., 7., 8.;
    const FreeVector freeVec(tempFreeVec);
    const Vector3D pos2 = position(freeVec);
    BOOST_CHECK_EQUAL(freeVec.head<3>(), pos2);
  }
}

template <size_t I>
struct functor {
  static constexpr size_t invoke() { return I * I * I; }
};

BOOST_AUTO_TEST_CASE(test_matrix_dimension_switch) {
  constexpr size_t imax = 20;
  for (size_t i = 0; i < imax; i++) {
    size_t val = template_switch<functor, 0, imax>(i);
    BOOST_CHECK_EQUAL(val, i * i * i);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
