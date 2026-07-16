// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cstddef>
#include <iomanip>
#include <memory>
#include <numbers>
#include <numeric>
#include <ostream>
#include <sstream>
#include <vector>

namespace Acts {

/// Class to implement pure virtual method of SurfaceBounds for testing only
class SurfaceBoundsStub : public SurfaceBounds {
 public:
  /// Implement ctor and pure virtual methods of SurfaceBounds
  explicit SurfaceBoundsStub(std::size_t nValues = 0) : m_values(nValues, 0) {
    std::iota(m_values.begin(), m_values.end(), 0);
  }

#if defined(__GNUC__) && (__GNUC__ == 13 || __GNUC__ == 14) && \
    !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  SurfaceBoundsStub(const SurfaceBoundsStub& other) = default;
  SurfaceBoundsStub& operator=(const SurfaceBoundsStub& other) = default;
#if defined(__GNUC__) && (__GNUC__ == 13 || __GNUC__ == 14) && \
    !defined(__clang__)
#pragma GCC diagnostic pop
#endif

  BoundsType type() const final { return eOther; }

  bool isCartesian() const final { return true; }

  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  std::vector<double> values() const final { return m_values; }

  bool inside(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return true;
  }

  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final {
    static_cast<void>(metric);
    return lposition;
  }

  Vector2 center() const final { return Vector2(0.0, 0.0); }

  bool inside(const Vector2& lposition,
              const BoundaryTolerance& boundaryTolerance) const final {
    static_cast<void>(lposition);
    static_cast<void>(boundaryTolerance);
    return true;
  }

  std::ostream& toStream(std::ostream& sl) const final {
    sl << "SurfaceBoundsStub";
    return sl;
  }

 private:
  std::vector<double> m_values;
};

}  // namespace Acts

using namespace Acts;

namespace ActsTests {

namespace {

struct StreamState {
  std::ios_base::fmtflags flags;
  std::streamsize precision;
  std::streamsize width;
  char fill;
};

StreamState setNonDefaultStreamState(std::ostringstream& stream) {
  stream << std::scientific << std::showpos << std::setfill('#')
         << std::setprecision(3);
  stream.width(17);
  return {stream.flags(), stream.precision(), stream.width(), stream.fill()};
}

void checkStreamState(const std::ostringstream& stream,
                      const StreamState& state) {
  BOOST_CHECK(stream.flags() == state.flags);
  BOOST_CHECK_EQUAL(stream.precision(), state.precision);
  BOOST_CHECK_EQUAL(stream.width(), state.width);
  BOOST_CHECK_EQUAL(stream.fill(), state.fill);
}

void checkStreamStatePreserved(const SurfaceBounds& bounds) {
  std::ostringstream stream;
  const auto state = setNonDefaultStreamState(stream);

  bounds.toStream(stream);

  checkStreamState(stream, state);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

/// Unit test for creating compliant/non-compliant SurfaceBounds object
BOOST_AUTO_TEST_CASE(SurfaceBoundsConstruction) {
  SurfaceBoundsStub u;
  SurfaceBoundsStub s(1);  // would act as std::size_t cast to SurfaceBounds
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

BOOST_AUTO_TEST_CASE(SurfaceBoundsToStreamPreservesStreamState) {
  std::vector<std::unique_ptr<SurfaceBounds>> bounds;
  bounds.push_back(
      std::make_unique<AnnulusBounds>(7.2, 12., 0.7, 1.3, Vector2{-2., 2.}));
  bounds.push_back(std::make_unique<ConeBounds>(std::numbers::pi / 8., 3., 6.));
  bounds.push_back(std::make_unique<CylinderBounds>(0.5, 10.));
  bounds.push_back(std::make_unique<DiamondBounds>(10., 20., 15., 5., 7.));
  bounds.push_back(std::make_unique<DiscTrapezoidBounds>(1., 5., 2., 6., 0.));
  bounds.push_back(std::make_unique<EllipseBounds>(1., 2., 3., 4.,
                                                   std::numbers::pi / 2., 0.));
  bounds.push_back(std::make_unique<LineBounds>(0.5, 20.));
  bounds.push_back(std::make_unique<RadialBounds>(1., 5.));
  bounds.push_back(std::make_unique<RectangleBounds>(10., 5.));
  bounds.push_back(std::make_unique<TrapezoidBounds>(1., 6., 2.));

  for (const auto& bound : bounds) {
    checkStreamStatePreserved(*bound);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
