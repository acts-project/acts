// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"
#include "Acts/Plugins/Detray/DetrayPayloadConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <memory>
#include <numbers>

#include <detray/io/frontend/payloads.hpp>

auto logger = Acts::getDefaultLogger("Test", Acts::Logging::INFO);

using namespace Acts;

BOOST_AUTO_TEST_SUITE(DetrayConversion)

BOOST_AUTO_TEST_CASE(DetrayTransformConversion) {
  auto transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(std::numbers::pi / 2., Vector3::UnitZ()));

  detray::io::transform_payload payload =
      DetrayPayloadConverter::convertTransform(transform);
  // Transform is correctly translated
  CHECK_CLOSE_ABS(payload.tr[0u], 1., std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.tr[1u], 2., std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.tr[2u], 3., std::numeric_limits<double>::epsilon());
  // Rotation is correctly translated
  RotationMatrix3 rotation = transform.rotation().transpose();
  CHECK_CLOSE_ABS(payload.rot[0u], rotation(0, 0),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[1u], rotation(0, 1),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[2u], rotation(0, 2),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[3u], rotation(1, 0),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[4u], rotation(1, 1),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[5u], rotation(1, 2),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[6u], rotation(2, 0),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[7u], rotation(2, 1),
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[8u], rotation(2, 2),
                  std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(DetrayMaskConversion) {
  // Test AnnulusBounds conversion
  {
    Acts::AnnulusBounds annulus(10., 20., 0.5, 1.5, {1., 2.}, 1.);
    auto payload = DetrayPayloadConverter::convertMask(annulus, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::annulus2);
    using enum detray::annulus2D::boundaries;
    CHECK_CLOSE_ABS(payload.boundaries[e_min_r], 10., 1e-10);  // min_r
    CHECK_CLOSE_ABS(payload.boundaries[e_max_r], 20., 1e-10);  // max_r
    CHECK_CLOSE_ABS(payload.boundaries[e_min_phi_rel], 0.5,
                    1e-10);  // min_phi_rel
    CHECK_CLOSE_ABS(payload.boundaries[e_max_phi_rel], 1.5,
                    1e-10);  // max_phi_rel
    CHECK_CLOSE_ABS(payload.boundaries[e_average_phi], 1.,
                    1e-10);                                     // average_phi
    CHECK_CLOSE_ABS(payload.boundaries[e_shift_x], 1., 1e-10);  // shift_x
    CHECK_CLOSE_ABS(payload.boundaries[e_shift_y], 2., 1e-10);  // shift_y
  }

  // Test RectangleBounds conversion
  {
    Acts::RectangleBounds rectangle(5., 10.);  // symmetric around origin
    auto payload = DetrayPayloadConverter::convertMask(rectangle, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::rectangle2);
    using enum detray::rectangle2D::boundaries;
    CHECK_CLOSE_ABS(payload.boundaries[e_half_x], 5., 1e-10);   // half_x
    CHECK_CLOSE_ABS(payload.boundaries[e_half_y], 10., 1e-10);  // half_y
  }

  // Test CylinderBounds conversion (normal and portal)
  {
    Acts::CylinderBounds cylinder(50., 100.);
    auto payload = DetrayPayloadConverter::convertMask(cylinder, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::cylinder2);
    using enum detray::concentric_cylinder2D::boundaries;
    CHECK_CLOSE_ABS(payload.boundaries[e_r], 50., 1e-10);          // r
    CHECK_CLOSE_ABS(payload.boundaries[e_lower_z], -100., 1e-10);  // lower_z
    CHECK_CLOSE_ABS(payload.boundaries[e_upper_z], 100., 1e-10);   // upper_z

    // Test portal case
    auto portalPayload = DetrayPayloadConverter::convertMask(cylinder, true);
    BOOST_CHECK(portalPayload.shape == detray::io::shape_id::portal_cylinder2);
    using enum detray::concentric_cylinder2D::boundaries;
    CHECK_CLOSE_ABS(portalPayload.boundaries[e_r], 50., 1e-10);  // r
    CHECK_CLOSE_ABS(portalPayload.boundaries[e_lower_z], -100.,
                    1e-10);  // lower_z
    CHECK_CLOSE_ABS(portalPayload.boundaries[e_upper_z], 100.,
                    1e-10);  // upper_z
  }

  // Test TrapezoidBounds conversion
  {
    Acts::TrapezoidBounds trapezoid(4., 8., 10.,
                                    0.);  // minHalfX, maxHalfX, halfY
    auto payload = DetrayPayloadConverter::convertMask(trapezoid, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::trapezoid2);
    using enum detray::trapezoid2D::boundaries;
    CHECK_CLOSE_ABS(payload.boundaries[e_half_length_0], 4.,
                    1e-10);  // half_length_0
    CHECK_CLOSE_ABS(payload.boundaries[e_half_length_1], 8.,
                    1e-10);  // half_length_1
    CHECK_CLOSE_ABS(payload.boundaries[e_half_length_2], 10.,
                    1e-10);  // half_length_2
    CHECK_CLOSE_ABS(payload.boundaries[e_divisor], 1. / (2. * 10.),
                    1e-10);  // divisor
  }

  // Test RadialBounds conversion
  {
    Acts::RadialBounds radial(5.,
                              15);  // minR, maxR, minPhi, maxPhi (full azimuth)
    auto payload = DetrayPayloadConverter::convertMask(radial, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::ring2);
    using enum detray::ring2D::boundaries;
    CHECK_CLOSE_ABS(payload.boundaries[e_inner_r], 5., 1e-10);   // inner_r
    CHECK_CLOSE_ABS(payload.boundaries[e_outer_r], 15., 1e-10);  // outer_r
  }
}

BOOST_AUTO_TEST_CASE(DetrayMaskConversionErrors) {
  // Test non-symmetric RectangleBounds throws
  // The regular rectangle bounds constructor does not allow this, but let's
  // test it anyway.
  {
    std::array<double, Acts::RectangleBounds::eSize> values = {
        -2., -3., 5., 4.};  // minX, minY, maxX, maxY
    Acts::RectangleBounds asymmetricRect(values);
    BOOST_CHECK_THROW(
        DetrayPayloadConverter::convertMask(asymmetricRect, false),
        std::runtime_error);
  }

  // Test partial azimuth RadialBounds throws
  {
    Acts::RadialBounds partialRadial(5., 15., 0.5 * std::numbers::pi);
    BOOST_CHECK_THROW(DetrayPayloadConverter::convertMask(partialRadial, false),
                      std::runtime_error);
  }

  {
    Acts::RadialBounds partialRadial(
        5., 15., std::numbers::pi,
        0.75 * std::numbers::pi);  // full azimuth, but internal rotation
    BOOST_CHECK_THROW(DetrayPayloadConverter::convertMask(partialRadial, false),
                      std::runtime_error);
  }

  // Test unsupported disc bounds type throws
  {
    class MockDiscBounds final : public Acts::SurfaceBounds {
     public:
      BoundsType type() const final { return BoundsType::eDisc; }
      std::vector<double> values() const final { return {}; }
      bool inside(const Acts::Vector2& /*loc*/,
                  const Acts::BoundaryTolerance& /*bcheck*/) const final {
        return false;
      }
      std::ostream& toStream(std::ostream& sl) const final { return sl; }
    };
    MockDiscBounds mockDisc;
    BOOST_CHECK_THROW(DetrayPayloadConverter::convertMask(mockDisc, false),
                      std::runtime_error);
  }

  // Test unknown bounds type returns unknown shape
  {
    class MockUnknownBounds final : public Acts::SurfaceBounds {
     public:
      BoundsType type() const final {
        return static_cast<BoundsType>(999);
      }  // Invalid type
      std::vector<double> values() const final { return {}; }
      bool inside(const Acts::Vector2& /*loc*/,
                  const Acts::BoundaryTolerance& /*bcheck*/) const final {
        return false;
      }
      std::ostream& toStream(std::ostream& sl) const final { return sl; }
    };
    MockUnknownBounds mockUnknown;
    auto payload = DetrayPayloadConverter::convertMask(mockUnknown, false);
    BOOST_CHECK(payload.shape == detray::io::shape_id::unknown);
  }
}

BOOST_AUTO_TEST_SUITE_END()
