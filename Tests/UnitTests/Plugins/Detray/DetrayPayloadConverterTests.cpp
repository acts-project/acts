// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <numbers>

#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/geometry_writer.hpp>
#include <detray/io/backend/homogeneous_material_reader.hpp>
#include <detray/io/backend/homogeneous_material_writer.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>
#include <detray/io/frontend/definitions.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/io/frontend/payloads.hpp>
#include <detray/io/json/json_io.hpp>
#include <detray/plugins/svgtools/illustrator.hpp>
#include <detray/plugins/svgtools/writer.hpp>
#include <detray/utils/consistency_checker.hpp>
#include <detray/utils/grid/detail/concepts.hpp>
#include <nlohmann/json_fwd.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/geometry/tracking_volume.hpp"

auto logger = Acts::getDefaultLogger("Test", Acts::Logging::INFO);

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetrayConversion)

BOOST_AUTO_TEST_CASE(DetrayTransformConversion) {
  auto transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(std::numbers::pi / 2., Vector3::UnitZ()));

  detray::io::transform_payload payload =
      DetrayConversionUtils::convertTransform(transform);
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
      bool isCartesian() const final { return true; }
      SquareMatrix2 boundToCartesianJacobian(
          const Vector2& /*lposition*/) const final {
        return SquareMatrix2::Identity();
      }
      SquareMatrix2 boundToCartesianMetric(
          const Vector2& /*lposition*/) const final {
        return SquareMatrix2::Identity();
      }
      std::vector<double> values() const final { return {}; }
      Vector2 center() const final { return Vector2::Zero(); }
      bool inside(const Vector2& /*lposition*/) const final { return false; }
      Vector2 closestPoint(const Vector2& lposition,
                           const SquareMatrix2& /*metric*/) const final {
        return lposition;
      }
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
      bool isCartesian() const final { return true; }
      SquareMatrix2 boundToCartesianJacobian(
          const Vector2& /*lposition*/) const final {
        return SquareMatrix2::Identity();
      }
      SquareMatrix2 boundToCartesianMetric(
          const Vector2& /*lposition*/) const final {
        return SquareMatrix2::Identity();
      }
      std::vector<double> values() const final { return {}; }
      Vector2 center() const final { return Vector2::Zero(); }
      bool inside(const Vector2& /*lposition*/) const final { return false; }
      Vector2 closestPoint(const Vector2& lposition,
                           const SquareMatrix2& /*metric*/) const final {
        return lposition;
      }
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

BOOST_AUTO_TEST_CASE(DetraySurfaceConversionTests) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  // Create a transform with translation and rotation
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(0.5 * std::numbers::pi, Vector3::UnitZ()));

  // Create rectangle bounds
  auto bounds = std::make_shared<RectangleBounds>(5., 10.);

  // Create a plane surface
  auto surface = Surface::makeShared<PlaneSurface>(transform, bounds);

  // Test surface conversion with Identifier strategy
  {
    DetrayPayloadConverter::Config cfg;
    cfg.sensitiveStrategy =
        DetrayPayloadConverter::Config::SensitiveStrategy::Identifier;

    DetrayPayloadConverter converter(cfg);

    // Test passive surface (default)
    {
      auto payload = converter.convertSurface(gctx, *surface);

      // Check type
      BOOST_CHECK(payload.type == detray::surface_id::e_passive);

      // Check transform
      CHECK_CLOSE_ABS(payload.transform.tr[0], 1., 1e-10);
      CHECK_CLOSE_ABS(payload.transform.tr[1], 2., 1e-10);
      CHECK_CLOSE_ABS(payload.transform.tr[2], 3., 1e-10);

      // Check mask
      BOOST_CHECK(payload.masks.at(0).shape ==
                  detray::io::shape_id::rectangle2);
      using enum detray::rectangle2D::boundaries;
      CHECK_CLOSE_ABS(payload.masks.at(0).boundaries[e_half_x], 5., 1e-10);
      CHECK_CLOSE_ABS(payload.masks.at(0).boundaries[e_half_y], 10., 1e-10);
    }

    // Test sensitive surface
    {
      auto sensitiveSurface =
          Surface::makeShared<PlaneSurface>(transform, bounds);
      sensitiveSurface->assignGeometryId(GeometryIdentifier().withSensitive(1));

      auto payload = converter.convertSurface(gctx, *sensitiveSurface);

      BOOST_CHECK(payload.type == detray::surface_id::e_sensitive);
    }
  }

  // Test surface conversion with DetectorElement strategy
  {
    DetrayPayloadConverter::Config cfg;
    cfg.sensitiveStrategy =
        DetrayPayloadConverter::Config::SensitiveStrategy::DetectorElement;

    DetrayPayloadConverter converter(cfg);

    // Test passive surface (no detector element)
    {
      auto payload = converter.convertSurface(gctx, *surface);
      BOOST_CHECK(payload.type == detray::surface_id::e_passive);
    }

    // Test sensitive surface with detector element
    {
      // Create detector element first
      auto detElement =
          std::make_shared<DetectorElementStub>(transform, bounds, 1.0);

      // Create surface using the detector element
      auto sensitiveSurface =
          Surface::makeShared<PlaneSurface>(bounds, *detElement);

      auto payload = converter.convertSurface(gctx, *sensitiveSurface);
      BOOST_CHECK(payload.type == detray::surface_id::e_sensitive);
    }
  }
}

BOOST_AUTO_TEST_CASE(DetrayPortalConversionTests) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  // Create a transform with translation and rotation
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(0.5 * std::numbers::pi, Vector3::UnitZ()));

  // Create rectangle bounds
  auto bounds = std::make_shared<RectangleBounds>(5., 10.);

  auto cvlBounds = std::make_shared<CylinderVolumeBounds>(5., 10., 10.);
  TrackingVolume tv(transform, cvlBounds);

  // Create a plane surface
  auto surface = Surface::makeShared<PlaneSurface>(transform, bounds);

  // Create a portal with the surface
  auto portal = std::make_shared<Portal>(Direction::AlongNormal(), surface, tv);

  // Test portal conversion
  {
    DetrayPayloadConverter::Config cfg;
    DetrayPayloadConverter converter(cfg);
    auto payload = converter.convertSurface(gctx, portal->surface(), true);

    // Portal should always be passive
    BOOST_CHECK(payload.type == detray::surface_id::e_portal);

    // Check transform is preserved
    CHECK_CLOSE_ABS(payload.transform.tr[0], 1., 1e-10);
    CHECK_CLOSE_ABS(payload.transform.tr[1], 2., 1e-10);
    CHECK_CLOSE_ABS(payload.transform.tr[2], 3., 1e-10);

    // Check mask - should be same as surface since rectangle doesn't have
    // special portal handling
    BOOST_CHECK(payload.masks.at(0).shape == detray::io::shape_id::rectangle2);
    using enum detray::rectangle2D::boundaries;
    CHECK_CLOSE_ABS(payload.masks.at(0).boundaries[e_half_x], 5., 1e-10);
    CHECK_CLOSE_ABS(payload.masks.at(0).boundaries[e_half_y], 10., 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(DetrayVolumeConversionTests) {
  // Create a transform with translation and rotation
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(0.5 * std::numbers::pi, Vector3::UnitZ()));

  DetrayPayloadConverter::Config cfg;
  DetrayPayloadConverter converter(cfg);

  // Test cylinder volume
  {
    auto cvlBounds = std::make_shared<CylinderVolumeBounds>(5., 10., 10.);
    auto volume =
        std::make_shared<TrackingVolume>(transform, cvlBounds, "TestCylinder");

    auto payload = converter.convertVolume(*volume);

    // Check type
    BOOST_CHECK(payload.type == detray::volume_id::e_cylinder);

    // Check name
    BOOST_CHECK_EQUAL(payload.name, "TestCylinder");

    // Check transform
    CHECK_CLOSE_ABS(payload.transform.tr[0], 1., 1e-10);
    CHECK_CLOSE_ABS(payload.transform.tr[1], 2., 1e-10);
    CHECK_CLOSE_ABS(payload.transform.tr[2], 3., 1e-10);
  }

  // Test cuboid volume
  {
    auto cuboidBounds = std::make_shared<CuboidVolumeBounds>(5., 10., 15.);
    auto volume =
        std::make_shared<TrackingVolume>(transform, cuboidBounds, "TestCuboid");

    auto payload = converter.convertVolume(*volume);

    BOOST_CHECK(payload.type == detray::volume_id::e_cuboid);
    BOOST_CHECK_EQUAL(payload.name, "TestCuboid");
  }

  // Test trapezoid volume
  {
    auto trapBounds = std::make_shared<TrapezoidVolumeBounds>(4., 8., 10., 15.);
    auto volume = std::make_shared<TrackingVolume>(transform, trapBounds,
                                                   "TestTrapezoid");

    auto payload = converter.convertVolume(*volume);

    BOOST_CHECK(payload.type == detray::volume_id::e_trapezoid);
    BOOST_CHECK_EQUAL(payload.name, "TestTrapezoid");
  }

  // Test unknown volume type
  {
    class MockVolumeBounds : public VolumeBounds {
     public:
      BoundsType type() const final { return BoundsType::eOther; }
      std::vector<double> values() const final { return {}; }
      bool inside(const Vector3& /*pos*/, double /*tol*/) const final {
        return false;
      }
      std::ostream& toStream(std::ostream& sl) const final { return sl; }

      std::vector<OrientedSurface> orientedSurfaces(
          const Transform3& /*trf*/) const final {
        return {};
      }

      Volume::BoundingBox boundingBox(const Transform3* /*trf*/,
                                      const Vector3& /*envelope*/,
                                      const Volume* /*entity*/) const final {
        return {nullptr, Vector3::Zero(), Vector3::Zero()};
      }
    };

    auto mockBounds = std::make_shared<MockVolumeBounds>();
    TrackingVolume tv(transform, mockBounds, "TestUnknown");

    auto volume =
        std::make_shared<TrackingVolume>(transform, mockBounds, "TestUnknown");

    auto payload = converter.convertVolume(*volume);

    BOOST_CHECK(payload.type == detray::volume_id::e_unknown);
    BOOST_CHECK_EQUAL(payload.name, "TestUnknown");
  }
}

// namespace {
//
// namespace detail {
//
// /// A functor that retrieves an acceleration struct and prints it
// struct accelerator_printer {
//   /// Print an acceleration structure
//   ///
//   /// @param accel_coll collection of acceleration structs
//   /// @param idx the specific grid to be checked
//   /// @param id type id of the material grid collection
//   template <typename accel_coll_t, typename index_t>
//   DETRAY_HOST void operator()(const accel_coll_t& accel_coll, const index_t
//   idx,
//                               std::stringstream& os) const {
//     // os << accel_coll[idx];
//   }
// };
//
// /// A functor that retrieves material and prints it
// struct material_printer {
//   /// Print material
//   ///
//   /// @param material_coll collection of material
//   /// @param idx the specific grid to be checked
//   /// @param id type id of the material grid collection
//   template <typename material_coll_t, typename index_t>
//   DETRAY_HOST void operator()(const material_coll_t& material_coll,
//                               const index_t idx, std::stringstream& os) const
//                               {
//     if constexpr (!detray::concepts::grid<
//                       typename material_coll_t::value_type>) {
//       os << material_coll[idx];
//     }
//   }
// };
//
// }  // namespace detail
//
// template <typename detector_t>
// inline std::string print_detector(
//     const detector_t& det, const typename detector_t::name_map& names = {}) {
//   // Gathers navigation information across navigator update calls
//   std::stringstream debug_stream{};
//
//   debug_stream << std::left << "[>] Detector " << det.name(names) << " has "
//                << det.volumes().size() << " volumes." << std::endl;
//
//   for (const auto [i, v_desc] : detray::views::enumerate(det.volumes())) {
//     detray::tracking_volume v{det, v_desc};
//
//     debug_stream << "[>>] Volume " << v.name(names) << std::endl;
//     debug_stream << v << std::endl;
//
//     debug_stream << "[>>>] Acceleration Structures:" << std::endl;
//     const auto acc_link = v_desc.accel_link();
//     for (std::size_t j = 0u; j < acc_link.size(); ++j) {
//       // An acceleration data structure link was set, but is invalid
//       if (!acc_link[j].is_invalid_id() && !acc_link[j].is_invalid_index()) {
//         debug_stream << j << ":" << std::endl;
//         det.accelerator_store().template visit<detail::accelerator_printer>(
//             acc_link[j], debug_stream);
//       }
//     }
//
//     debug_stream << "[>>>] Surfaces:" << std::endl;
//     for (const auto sf_desc : v.template surfaces<>()) {
//       detray::geometry::surface sf{det, sf_desc};
//       debug_stream << sf << std::endl;
//
//       // Check the surface material, if present
//       if (sf.has_material()) {
//         debug_stream << "[>>>>] Surface material:" << std::endl;
//         sf.template visit_material<detail::material_printer>(debug_stream);
//       }
//     }
//
//     debug_stream << std::endl;
//   }
//
//   return debug_stream.str();
// }
// }  // namespace

BOOST_AUTO_TEST_CASE(DetrayTrackingGeometryConversionTests) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto geoLogger = getDefaultLogger("Geo", Logging::VERBOSE);

  CylindricalTrackingGeometry cGeometry(gctx, true);
  auto tGeometry = cGeometry(*geoLogger);

  ObjVisualization3D obj;

  tGeometry->visualize(obj, gctx);

  obj.write("cylindrical.obj");

  vecmem::host_memory_resource mr;

  DetrayPayloadConverter::Config cfg;

  tGeometry->apply([&cfg](const TrackingVolume& volume) {
    if (volume.volumeName() == "Beampipe") {
      cfg.beampipeVolume = &volume;
    }
  });

  auto logger = getDefaultLogger("Cnv", Logging::DEBUG);
  DetrayPayloadConverter converter(cfg, std::move(logger));
  auto payloads = converter.convertTrackingGeometry(gctx, *tGeometry);

  const auto& detector = *payloads.detector;
  const auto& homogeneousMaterial = *payloads.homogeneousMaterial;
  auto& materialGrids = *payloads.materialGrids;
  const auto& surfaceGrids = *payloads.surfaceGrids;

  BOOST_CHECK_EQUAL(detector.volumes.size(), 6);
  for (const auto& volume : detector.volumes) {
    BOOST_CHECK(volume.type == detray::volume_id::e_cylinder);
  }
  BOOST_CHECK_EQUAL(detector.volumes.at(0).name, "Beampipe");
  BOOST_CHECK_EQUAL(detector.volumes.at(1).name, "World");
  BOOST_CHECK_EQUAL(detector.volumes.at(2).name, "L0");
  BOOST_CHECK_EQUAL(detector.volumes.at(3).name, "L1");
  BOOST_CHECK_EQUAL(detector.volumes.at(4).name, "L2");
  BOOST_CHECK_EQUAL(detector.volumes.at(5).name, "L3");

  // @HACK: At this time, the conversion introduces a number of dummy material slabs which
  //        should ultimately not be there.

  BOOST_CHECK_EQUAL(homogeneousMaterial.volumes.size(), 6);

  for (const auto& volume : detector.volumes) {
    auto it = std::ranges::find_if(
        homogeneousMaterial.volumes, [&](const auto& hMat) {
          return hMat.volume_link.link == volume.index.link;
        });

    if (it == homogeneousMaterial.volumes.end()) {
      BOOST_FAIL("No material slab found for volume: " + volume.name);
      continue;
    }

    auto& hMat = *it;

    // Currently, we expect exactly one slab per surface!
    BOOST_CHECK_EQUAL(volume.surfaces.size(), hMat.mat_slabs.size());
  }

  BOOST_CHECK_EQUAL(materialGrids.grids.size(), 2);

  BOOST_CHECK_EQUAL(surfaceGrids.grids.size(), 4);

  // Empirical binning config from construction
  const std::vector<std::pair<int, int>> kLayerBinning = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};

  for (const auto& [idx, it] : enumerate(surfaceGrids.grids)) {
    const auto [b0, b1] = kLayerBinning.at(idx);

    const auto [key, grids] = it;

    // We're only building one surface grid for each layer
    BOOST_CHECK_EQUAL(grids.size(), 1);

    auto& grid = grids.front();

    // All grids in this case are 2D
    BOOST_CHECK_EQUAL(grid.axes.size(), 2);

    BOOST_CHECK_EQUAL(b0, grid.axes.front().bins);
    BOOST_CHECK_EQUAL(b1, grid.axes.back().bins);

    std::size_t nBinsExp = grid.axes.front().bins * grid.axes.back().bins;
    BOOST_CHECK_EQUAL(grid.bins.size(), nBinsExp);
  }

  BOOST_CHECK_EQUAL(payloads.names.size(), 7);

  // Write payloads to JSON directly

  {
    detray::io::geo_header_payload header_data;
    header_data.common = detray::io::detail::basic_converter::to_payload(
        payloads.names.at(0), detray::io::geometry_writer::tag);
    header_data.sub_header.emplace();
    auto& geo_sub_header = header_data.sub_header.value();
    geo_sub_header.n_volumes = detector.volumes.size();
    geo_sub_header.n_surfaces = 0;
    for (const auto& volume : detector.volumes) {
      geo_sub_header.n_surfaces += volume.surfaces.size();
    }

    nlohmann::ordered_json out_json;
    out_json["header"] = header_data;
    out_json["data"] = detector;

    std::ofstream ofs{"Detector_geometry_direct.json"};
    ofs << out_json.dump(2) << std::endl;
  }

  {
    detray::io::homogeneous_material_header_payload header_data;
    header_data.common = detray::io::detail::basic_converter::to_payload(
        payloads.names.at(0), detray::io::homogeneous_material_writer::tag);
    header_data.sub_header.emplace();
    header_data.sub_header->n_rods = 0;
    header_data.sub_header->n_slabs = 0;

    for (const auto& hVol : homogeneousMaterial.volumes) {
      if (hVol.mat_rods.has_value()) {
        header_data.sub_header->n_rods += hVol.mat_rods->size();
      }
      header_data.sub_header->n_slabs += hVol.mat_slabs.size();
    }

    nlohmann::ordered_json out_json;
    out_json["header"] = header_data;
    out_json["data"] = homogeneousMaterial;

    std::ofstream ofs{"Detector_homogeneous_material_direct.json"};
    ofs << out_json.dump(2) << std::endl;
  }

  // Payloads DONE, let's actually build a detray detector from them.

  using detector_t =
      detray::detector<detray::default_metadata<detray::array<double>>>;

  // build detector
  detray::detector_builder<detector_t::metadata> detectorBuilder{};
  // (1) geometry
  detray::io::geometry_reader::from_payload<detector_t>(detectorBuilder,
                                                        detector);

  detray::io::homogeneous_material_reader::from_payload<detector_t>(
      detectorBuilder, homogeneousMaterial);

  detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
      from_payload<detector_t>(detectorBuilder, std::move(materialGrids));

  detray::io::surface_grid_reader<detector_t::surface_type,
                                  std::integral_constant<std::size_t, 0>,
                                  std::integral_constant<std::size_t, 2>>::
      from_payload<detector_t>(detectorBuilder, surfaceGrids);

  detector_t detrayDetector(detectorBuilder.build(mr));

  // Checks and print
  detray::detail::check_consistency(detrayDetector);

  detray::svgtools::illustrator illustrator(detrayDetector, payloads.names);
  illustrator.hide_eta_lines(true);
  illustrator.show_info(true);

  const auto svg_zr = illustrator.draw_detector(actsvg::views::z_r{});
  actsvg::style::stroke stroke_black = actsvg::style::stroke();
  auto zr_axis = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                        stroke_black, "z", "r");
  detray::svgtools::write_svg("test_svgtools_detector_zr", {
                                                               zr_axis,
                                                               svg_zr,
                                                           });

  const auto svg_xy = illustrator.draw_detector(actsvg::views::x_y{});
  auto xy_axis = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                        stroke_black, "x", "y");
  detray::svgtools::write_svg("test_svgtools_detector_xy", {
                                                               xy_axis,
                                                               svg_xy,
                                                           });

  auto writer_cfg = detray::io::detector_writer_config{}
                        .format(detray::io::format::json)
                        .replace_files(true);

  // std::cout << detray::utils::print_detector(detrayDetector, payloads.names)
  //           << std::endl;

  detray::io::write_detector(detrayDetector, payloads.names, writer_cfg);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
