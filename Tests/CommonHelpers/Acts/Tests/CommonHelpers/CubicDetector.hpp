// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"

#include <functional>
#include <vector>

namespace Acts {
namespace Test {

struct CubicDetector {
  /// Default constructor for the Cubic detector
  ///
  /// @param gctx the geometry context for this geometry at building time
  CubicDetector(const GeometryContext& gctx) : geoContext(gctx) {
    using namespace UnitLiterals;

    // Construct the rotation
    double rotationAngle = 90_degree;
    Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3 yPos(0., 1., 0.);
    Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Boundaries of the surfaces
    rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));

    // Material of the surfaces
    MaterialSlab matProp(makeBeryllium(), 0.5_mm);
    surfaceMaterial = std::make_shared<HomogeneousSurfaceMaterial>(matProp);
  }

  /// Call operator to build the standard cubic tracking geometry
  std::shared_ptr<const Acts::Experimental::Detector> operator()() {
    using namespace UnitLiterals;

    // Geometry Id generator
    Acts::Experimental::GeometryIdGenerator::Config generatorConfig;
    Acts::Experimental::GeometryIdGenerator generator(
        generatorConfig, Acts::getDefaultLogger("SequentialIdGenerator",
                                                Acts::Logging::VERBOSE));
    auto cache = generator.generateCache();

    // Set translation vectors
    double eps = 1_mm;
    std::vector<Vector3> translations;
    translations.push_back({-2_m, 0., 0.});
    translations.push_back({-1_m, 0., 0.});
    translations.push_back({1_m - eps, 0., 0.});
    translations.push_back({1_m + eps, 0., 0.});
    translations.push_back({2_m - eps, 0., 0.});
    translations.push_back({2_m + eps, 0., 0.});

    std::vector<double> rotAngle;
    rotAngle.push_back(0.);
    rotAngle.push_back(0.);
    rotAngle.push_back(0.026);
    rotAngle.push_back(-0.026);
    rotAngle.push_back(0.026);
    rotAngle.push_back(-0.026);

    // Construct surfaces
    std::array<std::shared_ptr<Surface>, 6> surfaces;
    for (unsigned int i = 0; i < translations.size(); i++) {
      RotationMatrix3 rotation_strip;
      double angle = rotAngle[i];
      Vector3 xPos(cos(angle), sin(angle), 0.);
      Vector3 yPos(-sin(angle), cos(angle), 0.);
      Vector3 zPos(0., 0., 1.);
      rotation_strip.col(0) = xPos;
      rotation_strip.col(1) = yPos;
      rotation_strip.col(2) = zPos;

      Transform3 trafo(Transform3::Identity() * rotation * rotation_strip);
      trafo.translation() = translations[i];

      // Create the detector element
      auto detElement = std::make_unique<DetectorElementStub>(
          trafo, rBounds, 1._um, surfaceMaterial);
      // And remember the surface
      surfaces[i] = detElement->surface().getSharedPtr();
      // Add it to the event store
      detectorStore.push_back(std::move(detElement));

      // Assign the geometry id
      generator.assignGeometryId(cache, *surfaces[i]);
    }

    // Build volume for surfaces with negative x-values
    Transform3 trafoVol1(Transform3::Identity());
    trafoVol1.translation() = Vector3(-1.5_m, 0., 0.);

    auto boundsVol1 = std::make_shared<CuboidVolumeBounds>(1.5_m, 0.5_m, 0.5_m);

    auto detectorVolume1 = Acts::Experimental::DetectorVolumeFactory::construct(
        Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
        "Volume1", trafoVol1, std::move(boundsVol1), {surfaces[0], surfaces[1]},
        std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
        Acts::Experimental::tryAllSubVolumes(),
        Acts::Experimental::tryAllPortalsAndSurfaces());

    // Build volume for surfaces with positive x-values
    Transform3 trafoVol2(Transform3::Identity());
    trafoVol2.translation() = Vector3(1.5_m, 0., 0.);

    std::vector<std::shared_ptr<Surface>> surfVec;
    for (unsigned int i = 2; i < 6; i++) {
      surfVec.push_back(surfaces[i]);
    }

    auto boundsVol2 = std::make_shared<CuboidVolumeBounds>(1.5_m, 0.5_m, 0.5_m);

    auto detectorVolume2 = Acts::Experimental::DetectorVolumeFactory::construct(
        Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
        "Volume2", trafoVol2, std::move(boundsVol2), surfVec,
        std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
        Acts::Experimental::tryAllSubVolumes(),
        Acts::Experimental::tryAllPortalsAndSurfaces());

    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        detectorVolumes = {detectorVolume1, detectorVolume2};

    generator.assignGeometryId(cache, *detectorVolume1);
    generator.assignGeometryId(cache, *detectorVolume2);

    // Connect the detector volumes
    auto portalContainer =
        Acts::Experimental::detail::CuboidalDetectorHelper::connect(
            geoContext, detectorVolumes, Acts::BinningValue::binX, {},
            Acts::Logging::VERBOSE);

    auto detector = Acts::Experimental::Detector::makeShared(
        "cubicDetector", detectorVolumes,
        Acts::Experimental::tryRootVolumes());  ///< Seems to be the default

    // Build and return detector
    return detector;
  }

  RotationMatrix3 rotation = RotationMatrix3::Identity();
  std::shared_ptr<const RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const ISurfaceMaterial> surfaceMaterial = nullptr;

  std::vector<std::unique_ptr<const DetectorElementStub>> detectorStore = {};

  std::reference_wrapper<const GeometryContext> geoContext;
};
}  // namespace Test
}  // namespace Acts
