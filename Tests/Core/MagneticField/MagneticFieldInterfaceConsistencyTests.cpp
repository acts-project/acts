// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file MagneticFieldInterfaceConsistencyTests.cpp
#define BOOST_TEST_MODULE Magnetic field interface consistency tests

// clang-format off
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
// clang-format on

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// This is the canonical interface that all field implementations
  /// need to comply with.
  /// This test group is should include ALL field implementations,
  /// to make sure they conform to the interface, even if they are
  /// not implicitly tested by some of the other tests (e.g. Propagation)
  /// The function does not assert any functionality, it just ensures
  /// the interface compiles
  template <class BField_t>
  void
  testInterfaceConsistency(const BField_t& field)
  {
    using Cache_t = typename BField_t::Cache;
    Vector3D pos(0, 0, 0);
    Vector3D B;
    ActsMatrixD<3, 3> gradient;

    // test interface method without cache
    field.getField(pos);
    field.getFieldGradient(pos, gradient);

    // test interface method with cache
    Cache_t cache;
    field.getField(pos, cache);
    field.getFieldGradient(pos, gradient, cache);
  }

  BOOST_AUTO_TEST_CASE(TestConstantBFieldInterfaceConsistency)
  {
    ConstantBField field(1, 1, 1);
    testInterfaceConsistency(field);
  }

  BOOST_AUTO_TEST_CASE(TestSolenoidBFieldInterfaceConsistency)
  {
    SolenoidBField field({100, 1000, 20, 5});
    testInterfaceConsistency(field);
  }

  BOOST_AUTO_TEST_CASE(TestInterpolatedBFieldMapInterfaceConsistency)
  {
    // define dummy mapper and field cell, we don't need them to do anything
    struct DummyFieldCell
    {
      Vector3D
      getField(const Vector3D&) const
      {
        return {0, 0, 0};
      }
      bool
      isInside(const Vector3D&) const
      {
        return true;
      }
    };

    struct DummyMapper : DummyFieldCell
    {
      concept::AnyFieldCell<>
      getFieldCell(const Vector3D&) const
      {
        return DummyFieldCell();
      }
      std::vector<size_t>
      getNBins() const
      {
        return {42};
      }
      std::vector<double>
      getMin() const
      {
        return {5};
      }
      std::vector<double>
      getMax() const
      {
        return {15};
      }
    };

    InterpolatedBFieldMap::Config config;
    config.scale  = 1.;
    config.mapper = DummyMapper();

    // create BField service
    InterpolatedBFieldMap b(std::move(config));

    testInterfaceConsistency(b);
  }

  BOOST_AUTO_TEST_CASE(TestSharedBFieldInterfaceConsistency)
  {
    SharedBField<ConstantBField> field(
        std::make_shared<ConstantBField>(Vector3D(1, 1, 1)));
    testInterfaceConsistency(field);
  }
}  // namespace Test

}  // namespace Acts
