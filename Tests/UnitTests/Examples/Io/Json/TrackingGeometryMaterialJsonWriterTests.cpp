// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "ActsExamples/Io/Json/TrackingGeometryMaterialJsonWriter.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <numbers>

using namespace Acts;
using namespace ActsExamples;

BOOST_AUTO_TEST_CASE(ExportLegacyProtoMaterialPlaceholderRanges) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  world->assignGeometryId(GeometryIdentifier().withVolume(1));
  auto surface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  // DD4hep uses zero-width ranges until the surface bounds are available.
  BinUtility bins(8, -std::numbers::pi, std::numbers::pi, closed,
                  AxisDirection::AxisPhi);
  bins += BinUtility(3, 0., 0., open, AxisDirection::AxisZ);
  auto proto = std::make_shared<ProtoSurfaceMaterial>(
      bins, MappingType::Default, "cylinder");
  surface->assignSurfaceMaterial(proto);
  world->addSurface(surface);
  TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);

  ActsTests::TemporaryDirectory tmp;
  for (const auto* extension : {".json", ".cbor"}) {
    TrackingGeometryMaterialJsonWriter::Config config;
    config.filePath = tmp.path() / (std::string("material") + extension);
    TrackingGeometryMaterialJsonWriter(config, Logging::WARNING)
        .write(GeometryContext::dangerouslyDefaultConstruct(), geometry);
    const auto material =
        TrackingGeometryMaterialJsonConverter{}.fromFile(config.filePath);
    const auto* decoded = dynamic_cast<const ProtoSurfaceMaterial*>(
        material.keyedSurfaces.at("cylinder").material.get());
    BOOST_REQUIRE(decoded != nullptr);
    const auto& axes = decoded->binning().binningData();
    BOOST_REQUIRE_EQUAL(axes.size(), 2);
    BOOST_CHECK_EQUAL(axes[0].bins(), 8);
    BOOST_CHECK_EQUAL(axes[1].bins(), 3);
    BOOST_CHECK_EQUAL(axes[1].min, -100.);
    BOOST_CHECK_EQUAL(axes[1].max, 100.);
  }
  // Export must leave the geometry's original placeholder untouched.
  BOOST_CHECK(surface->surfaceMaterialSharedPtr() == proto);
  BOOST_CHECK_EQUAL(proto->binning().binningData()[1].min, 0.);
  BOOST_CHECK_EQUAL(proto->binning().binningData()[1].max, 0.);
}
