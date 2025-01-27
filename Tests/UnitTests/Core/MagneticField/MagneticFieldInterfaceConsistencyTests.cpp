// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"

namespace Acts {
namespace Test {

// Create a test context
MagneticFieldContext mfContext = MagneticFieldContext();

/// This is the canonical interface that all field implementations
/// need to comply with.
/// This test group is should include ALL field implementations,
/// to make sure they conform to the interface, even if they are
/// not implicitly tested by some of the other tests (e.g. Propagation)
/// The function does not assert any functionality, it just ensures
/// the interface compiles
template <class BField_t>
void testInterfaceConsistency(const BField_t& field) {
  using Cache_t = typename BField_t::Cache;
  Vector3 pos(0, 0, 0);
  Vector3 B;
  ActsMatrix<3, 3> gradient;

  // test interface method without cache
  field.getField(pos);
  field.getFieldGradient(pos, gradient);

  // test interface method with cache
  Cache_t cache(mfContext);
  field.getField(pos, cache);
  field.getFieldGradient(pos, gradient, cache);
}

BOOST_AUTO_TEST_CASE(TestConstantBFieldInterfaceConsistency) {
  ConstantBField field(1, 1, 1);
  testInterfaceConsistency(field);
}

BOOST_AUTO_TEST_CASE(TestSolenoidBFieldInterfaceConsistency) {
  SolenoidBField field({100, 1000, 20, 5});
  testInterfaceConsistency(field);
}

BOOST_AUTO_TEST_CASE(TestInterpolatedBFieldMapInterfaceConsistency) {
  // define dummy mapper and field cell, we don't need them to do anything
  struct DummyFieldCell {
    Vector3 getField(const Vector3&) const { return {0, 0, 0}; }
    bool isInside(const Vector3&) const { return true; }
  };

  struct DummyMapper : DummyFieldCell {
    using FieldCell = DummyFieldCell;

    DummyFieldCell getFieldCell(const Vector3&) const {
      return DummyFieldCell();
    }
    std::vector<std::size_t> getNBins() const { return {42}; }
    std::vector<double> getMin() const { return {5}; }
    std::vector<double> getMax() const { return {15}; }
  };

  DummyMapper m;
  InterpolatedBFieldMap<DummyMapper>::Config config(std::move(m));
  config.scale = 1.;

  // create BField service
  InterpolatedBFieldMap<DummyMapper> b(std::move(config));

  testInterfaceConsistency(b);
}

}  // namespace Test
}  // namespace Acts
