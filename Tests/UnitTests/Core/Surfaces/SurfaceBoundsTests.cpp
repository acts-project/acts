// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <cstddef>
#include <numeric>
#include <ostream>
#include <vector>

namespace Acts {

/// Class to implement pure virtual method of SurfaceBounds for testing only
class SurfaceBoundsStub : public SurfaceBounds {
 public:
  /// Implement ctor and pure virtual methods of SurfaceBounds
  explicit SurfaceBoundsStub(size_t nValues = 0) : m_values(nValues, 0) {
    std::iota(m_values.begin(), m_values.end(), 0);
  }

#if defined(__GNUC__) && __GNUC__ == 13 && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  SurfaceBoundsStub(const SurfaceBoundsStub& other) = default;
#if defined(__GNUC__) && __GNUC__ == 13 && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

  ~SurfaceBoundsStub() override = default;
  BoundsType type() const final {
    return SurfaceBounds::eOther;
  }
  std::vector<double> values() const override {
    return m_values;
  }
  bool inside(const Vector2& /*lpos*/,
              const BoundaryCheck& /*bcheck*/) const final {
    return true;
  }

  std::ostream& toStream(std::ostream& sl) const final {
    sl << "SurfaceBoundsStub";
    return sl;
  }

 private:
  std::vector<double> m_values;
};

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant SurfaceBounds object
BOOST_AUTO_TEST_CASE(SurfaceBoundsConstruction) {
  SurfaceBoundsStub u;
  SurfaceBoundsStub s(1);  // would act as size_t cast to SurfaceBounds
  SurfaceBoundsStub t(s);
  SurfaceBoundsStub v(u);
}
BOOST_AUTO_TEST_CASE(SurfaceBoundsProperties) {
  SurfaceBoundsStub surface(5);
  std::vector<double> reference{0, 1, 2, 3, 4};
  const auto& boundValues = surface.values();
  BOOST_CHECK_EQUAL_COLLECTIONS(reference.cbegin(), reference.cend(),
                                boundValues.cbegin(), boundValues.cend());
}
/// Unit test for testing SurfaceBounds properties
BOOST_AUTO_TEST_CASE(SurfaceBoundsEquality) {
  SurfaceBoundsStub surface(1);
  SurfaceBoundsStub copiedSurface(surface);
  SurfaceBoundsStub differentSurface(2);
  BOOST_CHECK_EQUAL(surface, copiedSurface);
  BOOST_CHECK_NE(surface, differentSurface);
  SurfaceBoundsStub assignedSurface;
  assignedSurface = surface;
  BOOST_CHECK_EQUAL(surface, assignedSurface);
  const auto& surfaceboundValues = surface.values();
  const auto& assignedboundValues = assignedSurface.values();
  BOOST_CHECK_EQUAL_COLLECTIONS(
      surfaceboundValues.cbegin(), surfaceboundValues.cend(),
      assignedboundValues.cbegin(), assignedboundValues.cend());
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
