// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;

namespace ActssTests {

// Create a test context
MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_SUITE(MagneticFieldSuite)

BOOST_AUTO_TEST_CASE(TypeErasedCacheType) {
  bool constructor_called = false;
  bool destructor_called = false;

  struct MyCache {
    MyCache(int value, bool* ctor, bool* dtor) : m_value{value}, m_dtor{dtor} {
      (*ctor) = true;
    }
    ~MyCache() { (*m_dtor) = true; }
    int m_value;
    bool* m_dtor;
  };

  BOOST_CHECK(!constructor_called);
  BOOST_CHECK(!destructor_called);

  {
    MagneticFieldProvider::Cache cache{
        MagneticFieldProvider::Cache(std::in_place_type<MyCache>, 42,
                                     &constructor_called, &destructor_called)};
    BOOST_CHECK(constructor_called);
    BOOST_CHECK(!destructor_called);

    MyCache& v = cache.as<MyCache>();
    BOOST_CHECK_EQUAL(v.m_value, 42);
    v.m_value = 65;

    MyCache& v2 = cache.as<MyCache>();
    BOOST_CHECK_EQUAL(v2.m_value, 65);
  }

  BOOST_CHECK(constructor_called);
  BOOST_CHECK(destructor_called);
}

/// Field that changes linearly with the position, B = B0 + G x. It does not
/// implement the gradient, so it uses the default of the base class.
class LinearField final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  LinearField(Vector3 b0, SquareMatrix3 g) : m_b0(std::move(b0)), m_g(g) {}

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& position,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    return Result<Vector3>::success(m_b0 + m_g * position);
  }

 private:
  Vector3 m_b0;
  SquareMatrix3 m_g;
};

BOOST_AUTO_TEST_CASE(FieldGradientDefault) {
  SquareMatrix3 g;
  g << 1., 2., 3., 4., 5., 6., 7., 8., 9.;
  const LinearField field(Vector3(0.1, -0.2, 2.), 1e-3 * g);
  auto cache = field.makeCache(mfContext);

  BOOST_CHECK(!field.providesFieldGradient());
  auto res = field.getFieldAndGradient(Vector3(1., 2., 3.), cache);
  BOOST_REQUIRE(!res.ok());
  BOOST_CHECK(res.error() == MagneticFieldError::NotImplemented);
}

BOOST_AUTO_TEST_CASE(FieldGradientNumerically) {
  SquareMatrix3 g;
  g << 1., 2., 3., 4., 5., 6., 7., 8., 9.;
  g *= 1e-3;
  const Vector3 b0(0.1, -0.2, 2.);
  const LinearField field(b0, g);
  auto cache = field.makeCache(mfContext);

  const Vector3 position(10., -20., 30.);
  auto res = getFieldAndGradientNumerically(field, position, cache, 1.);
  BOOST_REQUIRE(res.ok());
  CHECK_CLOSE_ABS(res->field, Vector3(b0 + g * position), 1e-15);
  // Central differences are exact for a linear field up to round-off
  CHECK_CLOSE_ABS(res->gradient, g, 1e-14);
}

BOOST_AUTO_TEST_CASE(FieldGradientConstant) {
  const Vector3 b(0.1, -0.2, 2.);
  const Vector3 position(10., -20., 30.);

  const ConstantBField constant(b);
  auto constantCache = constant.makeCache(mfContext);
  BOOST_CHECK(constant.providesFieldGradient());
  auto constantRes = constant.getFieldAndGradient(position, constantCache);
  BOOST_REQUIRE(constantRes.ok());
  BOOST_CHECK_EQUAL(constantRes->field, b);
  BOOST_CHECK_EQUAL(constantRes->gradient, SquareMatrix3::Zero());

  const NullBField null;
  auto nullCache = null.makeCache(mfContext);
  BOOST_CHECK(null.providesFieldGradient());
  auto nullRes = null.getFieldAndGradient(position, nullCache);
  BOOST_REQUIRE(nullRes.ok());
  BOOST_CHECK_EQUAL(nullRes->field, Vector3::Zero());
  BOOST_CHECK_EQUAL(nullRes->gradient, SquareMatrix3::Zero());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActssTests
