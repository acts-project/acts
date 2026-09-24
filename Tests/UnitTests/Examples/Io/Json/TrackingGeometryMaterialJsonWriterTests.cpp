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

using namespace Acts;
using namespace ActsExamples;

BOOST_AUTO_TEST_CASE(ExportPreservesDeferredProtoMaterialRanges) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  world->assignGeometryId(GeometryIdentifier().withVolume(1));
  auto surface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  auto proto = std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D(
          {AxisSpec::DeferredEquidistant(8, AxisDirection::AxisRPhi),
           AxisSpec::DeferredEquidistant(3, AxisDirection::AxisZ)}),
      MappingType::Default, "cylinder");
  surface->assignSurfaceMaterial(proto);
  world->addSurface(surface);
  TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);

  ActsTests::TemporaryDirectory tmp;
  for (const auto* extension : {".json", ".cbor"}) {
    TrackingGeometryMaterialJsonWriter::Config config;
    config.filePath = tmp.path() / (std::string("material") + extension);
    TrackingGeometryMaterialJsonWriter(config, Logging::WARNING)
        .write(geometry);
    const auto material =
        TrackingGeometryMaterialJsonConverter{}.fromFile(config.filePath);
    const auto* decoded = dynamic_cast<const ProtoGridSurfaceMaterial*>(
        material.keyedSurfaces.at("cylinder").material.get());
    BOOST_REQUIRE(decoded != nullptr);
    const auto& axes = decoded->binning().axisSpecs();
    BOOST_REQUIRE_EQUAL(axes.size(), 2);
    BOOST_CHECK_EQUAL(axes[0].asEquidistant().nBins, 8);
    BOOST_CHECK_EQUAL(axes[1].asEquidistant().nBins, 3);
    for (const auto& axis : axes) {
      BOOST_CHECK(!axis.asEquidistant().min);
      BOOST_CHECK(!axis.asEquidistant().max);
    }
  }
  // Export must leave the geometry's original placeholder untouched.
  BOOST_CHECK(surface->surfaceMaterialSharedPtr() == proto);
  BOOST_CHECK(proto->binning().isDeferred());
}
