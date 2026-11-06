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
#include "Acts/MagneticField/MagneticFieldQuery.hpp"
#include "Acts/MagneticField/NullBField.hpp"

#include <memory>
#include <type_traits>
#include <utility>

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

/// Field that does not implement the gradient, so it uses the default of the
/// base class
class FieldWithoutGradient final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& /*position*/,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    return Result<Vector3>::success(Vector3(0., 0., 2.));
  }
};

BOOST_AUTO_TEST_CASE(FieldGradientDefault) {
  const FieldWithoutGradient field;
  auto cache = field.makeCache(mfContext);

  BOOST_CHECK(!field.providesFieldGradient());
  auto res = field.getFieldAndGradient(Vector3(1., 2., 3.), cache);
  BOOST_REQUIRE(!res.ok());
  BOOST_CHECK(res.error() == MagneticFieldError::NotImplemented);
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

namespace {
class CountingField final : public Acts::MagneticFieldProvider {
 public:
  struct LookupCache {
    const CountingField* owner;
    int lookups = 0;
  };
  mutable int cacheCreations = 0;
  mutable const Acts::MagneticFieldContext* lastContext = nullptr;

  Cache makeCache(const Acts::MagneticFieldContext& mctx) const override {
    lastContext = &mctx;
    ++cacheCreations;
    return Cache(std::in_place_type<LookupCache>, this, 0);
  }

  Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                       Cache& cache) const override {
    auto& state = cache.as<LookupCache>();
    BOOST_CHECK(state.owner == this);
    if (position.x() < 0.) {
      return Acts::Result<Acts::Vector3>::failure(
          Acts::MagneticFieldError::OutOfBounds);
    }
    return Acts::Result<Acts::Vector3>::success(
        Acts::Vector3(0., 0., ++state.lookups));
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(FieldQueryOwnsAndReusesCache) {
  CountingField provider;
  Acts::MagneticFieldContext context;
  Acts::MagneticFieldQuery first(provider, context);
  BOOST_CHECK_EQUAL(provider.cacheCreations, 1);
  for (int i = 1; i <= 3; ++i) {
    auto result = first.getField(Acts::Vector3::Zero());
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK_EQUAL(result->z(), i);
  }
  BOOST_CHECK_EQUAL(provider.cacheCreations, 1);

  // A second operation has independent cache state, even with one provider.
  Acts::MagneticFieldQuery second(provider, context);
  BOOST_CHECK_EQUAL(provider.cacheCreations, 2);
  auto result = second.getField(Acts::Vector3::Zero());
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->z(), 1.);
  auto error = first.getField(Acts::Vector3(-1., 0., 0.));
  BOOST_REQUIRE(!error.ok());
  BOOST_CHECK(error.error() == Acts::MagneticFieldError::OutOfBounds);
}

// Borrowing an lvalue is supported; temporaries must not silently dangle.
using FieldProviderRef = const Acts::MagneticFieldProvider&;
static_assert(
    std::is_constructible_v<Acts::MagneticFieldQuery, FieldProviderRef,
                            const Acts::MagneticFieldContext&>);
static_assert(
    !std::is_constructible_v<Acts::MagneticFieldQuery, FieldProviderRef,
                             Acts::MagneticFieldContext&&>);
static_assert(
    !std::is_constructible_v<Acts::MagneticFieldQuery, FieldProviderRef,
                             const Acts::MagneticFieldContext&&>);

static_assert(
    !std::is_constructible_v<Acts::MagneticFieldQuery, CountingField&&,
                             const Acts::MagneticFieldContext&>);
static_assert(
    !std::is_constructible_v<Acts::MagneticFieldQuery, const CountingField&&,
                             const Acts::MagneticFieldContext&>);

BOOST_AUTO_TEST_CASE(FieldQueryRetainsContextAcrossMoves) {
  CountingField provider;
  const Acts::MagneticFieldContext firstContext;
  const Acts::MagneticFieldContext secondContext;
  Acts::MagneticFieldQuery first(provider, firstContext);
  BOOST_CHECK(provider.lastContext == &first.context());
  BOOST_CHECK(&first.context() == &firstContext);
  BOOST_REQUIRE(first.getField(Acts::Vector3::Zero()).ok());

  Acts::MagneticFieldQuery moved(std::move(first));
  BOOST_CHECK(&moved.context() == &firstContext);
  BOOST_CHECK_EQUAL(provider.cacheCreations, 1);

  CountingField secondProvider;
  Acts::MagneticFieldQuery second(secondProvider, secondContext);
  BOOST_CHECK(secondProvider.lastContext == &second.context());
  BOOST_CHECK(&second.context() == &secondContext);
  second = std::move(moved);
  BOOST_CHECK(&second.context() == &firstContext);
  BOOST_CHECK_EQUAL(provider.cacheCreations, 1);
  BOOST_CHECK_EQUAL(secondProvider.cacheCreations, 1);
  auto result = second.getField(Acts::Vector3::Zero());
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->z(), 2.);
}
