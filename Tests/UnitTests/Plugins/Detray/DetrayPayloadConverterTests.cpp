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
#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"
#include "Acts/Plugins/Detray/DetrayPayloadConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <numbers>

#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>
#include <detray/io/frontend/definitions.hpp>
#include <detray/io/frontend/payloads.hpp>
#include <detray/plugins/svgtools/illustrator.hpp>
#include <detray/plugins/svgtools/writer.hpp>
#include <detray/utils/consistency_checker.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

auto logger = Acts::getDefaultLogger("Test", Acts::Logging::INFO);

using namespace Acts;
using namespace Acts::Test;

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

BOOST_AUTO_TEST_CASE(DetraySurfaceConversionTests) {
  GeometryContext gctx;

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
      BOOST_CHECK(payload.mask.shape == detray::io::shape_id::rectangle2);
      using enum detray::rectangle2D::boundaries;
      CHECK_CLOSE_ABS(payload.mask.boundaries[e_half_x], 5., 1e-10);
      CHECK_CLOSE_ABS(payload.mask.boundaries[e_half_y], 10., 1e-10);
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
      auto detElement = std::make_shared<Acts::Test::DetectorElementStub>(
          transform, bounds, 1.0);

      // Create surface using the detector element
      auto sensitiveSurface =
          Surface::makeShared<PlaneSurface>(bounds, *detElement);

      auto payload = converter.convertSurface(gctx, *sensitiveSurface);
      BOOST_CHECK(payload.type == detray::surface_id::e_sensitive);
    }
  }
}

BOOST_AUTO_TEST_CASE(DetrayPortalConversionTests) {
  GeometryContext gctx;

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
    BOOST_CHECK(payload.mask.shape == detray::io::shape_id::rectangle2);
    using enum detray::rectangle2D::boundaries;
    CHECK_CLOSE_ABS(payload.mask.boundaries[e_half_x], 5., 1e-10);
    CHECK_CLOSE_ABS(payload.mask.boundaries[e_half_y], 10., 1e-10);
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

BOOST_AUTO_TEST_CASE(DetrayTrackingGeometryConversionTests) {
  GeometryContext gctx;

  CylindricalTrackingGeometry cGeometry(gctx, true);
  auto tGeometry = cGeometry();

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

  using detector_t =
      detray::detector<detray::default_metadata<detray::array<double>>>;

  detector_t::name_map names = {{0u, "Detector"}};

  // build detector
  detray::detector_builder<detector_t::metadata> detectorBuilder{};
  // (1) geometry
  detray::io::geometry_reader::from_payload<detector_t>(detectorBuilder, names,
                                                        *payloads.detector);

  detray::io::homogeneous_material_reader::from_payload<detector_t>(
      detectorBuilder, names, *payloads.homogeneousMaterial);

  detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
      from_payload<detector_t>(detectorBuilder, names,
                               std::move(*payloads.materialGrids));

#if 0
  // (2b) material grids
  if constexpr (detray::concepts::has_material_maps<detector_t>) {
    if (options.convertMaterial) {
      detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                         detray::io::material_id>
          materialGridsPayload =
              DetrayMaterialConverter::convertGridSurfaceMaterial(
                  cCache, detector, logger());
      detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
          from_payload<detector_t>(detectorBuilder, names,
                                   std::move(materialGridsPayload));
    }
  }

  // (3) surface grids
  if (options.convertSurfaceGrids) {
    detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>
        surfaceGridsPayload =
            DetraySurfaceGridsConverter::convertSurfaceGrids(detector);

    // Capacity 0 (dynamic bin size) and dimension 2 (2D grids)
    detray::io::surface_grid_reader<typename detector_t::surface_type,
                                    std::integral_constant<std::size_t, 0>,
                                    std::integral_constant<std::size_t, 2>>::
        template from_payload<detector_t>(detectorBuilder, names,
                                          std::move(surfaceGridsPayload));
  }

#endif

  detector_t detrayDetector(detectorBuilder.build(mr));

  // Checks and print
  detray::detail::check_consistency(detrayDetector);

  detray::svgtools::illustrator illustrator(detrayDetector, names);
  illustrator.hide_eta_lines(true);
  illustrator.show_info(true);

  illustrator.draw_detector(actsvg::views::z_r{});
  const auto svg_zr = illustrator.draw_detector(actsvg::views::z_r{});
  actsvg::style::stroke stroke_black = actsvg::style::stroke();
  auto zr_axis = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                        stroke_black, "z", "r");
  detray::svgtools::write_svg("test_svgtools_detector_zr", {
                                                               zr_axis,
                                                               svg_zr,
                                                           });

  auto writer_cfg = detray::io::detector_writer_config{}
                        .format(detray::io::format::json)
                        .replace_files(true);

  detray::io::write_detector(detrayDetector, names, writer_cfg);

  // // Check payload size
  // BOOST_CHECK_EQUAL(payload.volumes.size(), 1);
  // BOOST_CHECK_EQUAL(payload.volumes[0].type, detray::volume_id::e_cylinder);
  // BOOST_CHECK_EQUAL(payload.volumes[0].name, "CylinderVolume");
}

BOOST_AUTO_TEST_SUITE_END()
