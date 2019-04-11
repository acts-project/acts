// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <functional>
#include <vector>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/detail/DefaultDetectorElementBase.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

struct CubicTrackingGeometry {
  /// Default constructor for the Cubit tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  CubicTrackingGeometry(std::reference_wrapper<const GeometryContext> gctx)
      : geoContext(gctx) {
    using namespace UnitLiterals;

    // Construct the rotation
    double rotationAngle = 90_degree;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Boundaries of the surfaces
    rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));

    // Material of the surfaces
    MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5_mm);
    surfaceMaterial = std::shared_ptr<const ISurfaceMaterial>(
        new HomogeneousSurfaceMaterial(matProp));
  }
<<<<<<< HEAD

  /// Call operator to build the standard cubic tracking geometry
  std::shared_ptr<const TrackingGeometry> operator()() {
    using namespace UnitLiterals;

    // Set translation vectors
    double eps = 1_mm;
    std::vector<Vector3D> translations;
    translations.push_back({-2_m, 0., 0.});
    translations.push_back({-1_m, 0., 0.});
    translations.push_back({1_m - eps, 0., 0.});
    translations.push_back({1_m + eps, 0., 0.});
    translations.push_back({2_m - eps, 0., 0.});
    translations.push_back({2_m + eps, 0., 0.});

    // Construct surfaces
    std::array<std::shared_ptr<const Surface>, 6> surfaces;
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];
      // Create the detector element
      auto detElement = std::make_unique<const DetectorElementStub>(
          std::make_shared<const Transform3D>(trafo), rBounds, 1._um,
          surfaceMaterial);
      // And remember the surface
      surfaces[i] = detElement->surface().getSharedPtr();
      // Add it to the event store
      detectorStore.push_back(std::move(detElement));
    }

    // Construct layers
    std::array<LayerPtr, 6> layers;
    for (i = 0; i < 6; i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];

      std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfaces[i]));

      layers[i] = PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
                                     rBounds, std::move(surArray), 1._mm);

      auto mutableSurface = const_cast<Surface*>(surfaces[i].get());
      mutableSurface->associateLayer(*layers[i]);
    }

    // Build volume for surfaces with negative x-values
    Transform3D trafoVol1(Transform3D::Identity());
    trafoVol1.translation() = Vector3D(-1.5_m, 0., 0.);

    auto boundsVol =
        std::make_shared<const CuboidVolumeBounds>(1.5_m, 0.5_m, 0.5_m);

    LayerArrayCreator::Config lacConfig;
    LayerArrayCreator layArrCreator(
        lacConfig, getDefaultLogger("LayerArrayCreator", Logging::INFO));

    LayerVector layVec;
    layVec.push_back(layers[0]);
    layVec.push_back(layers[1]);
    std::unique_ptr<const LayerArray> layArr1(layArrCreator.layerArray(
        geoContext, layVec, -2_m - 1._mm, -1._m + 1._mm, BinningType::arbitrary,
        BinningValue::binX));

    auto trackVolume1 = TrackingVolume::create(
        std::make_shared<const Transform3D>(trafoVol1), boundsVol, nullptr,
        std::move(layArr1), nullptr, {}, "Volume 1");
    trackVolume1->sign(GeometrySignature::Global);

    // Build volume for surfaces with positive x-values
    Transform3D trafoVol2(Transform3D::Identity());
    trafoVol2.translation() = Vector3D(1.5_m, 0., 0.);

    layVec.clear();
    for (i = 2; i < 6; i++)
      layVec.push_back(layers[i]);
    std::unique_ptr<const LayerArray> layArr2(
        layArrCreator.layerArray(geoContext, layVec, 1._m - 2._mm, 2._m + 2._mm,
                                 BinningType::arbitrary, BinningValue::binX));

    auto trackVolume2 = TrackingVolume::create(
        std::make_shared<const Transform3D>(trafoVol2), boundsVol, nullptr,
        std::move(layArr2), nullptr, {}, "Volume 2");

    trackVolume2->sign(GeometrySignature::Global);

    // Glue volumes
    trackVolume2->glueTrackingVolume(
        geoContext, BoundarySurfaceFace::negativeFaceYZ, trackVolume1.get(),
        BoundarySurfaceFace::positiveFaceYZ);

    trackVolume1->glueTrackingVolume(
        geoContext, BoundarySurfaceFace::positiveFaceYZ, trackVolume2.get(),
        BoundarySurfaceFace::negativeFaceYZ);

    // Build world volume
    Transform3D trafoWorld(Transform3D::Identity());
    trafoWorld.translation() = Vector3D(0., 0., 0.);

    auto worldVol =
        std::make_shared<const CuboidVolumeBounds>(3._m, 0.5_m, 0.5_m);

    std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

    tapVec.push_back(std::make_pair(trackVolume1, Vector3D(-1.5_m, 0., 0.)));
    tapVec.push_back(std::make_pair(trackVolume2, Vector3D(1.5_m, 0., 0.)));

    std::vector<float> binBoundaries = {-3._m, 0., 3._m};

    BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
    std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

    std::shared_ptr<const TrackingVolumeArray> trVolArr(
        new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

    MutableTrackingVolumePtr mtvpWorld(
        TrackingVolume::create(std::make_shared<const Transform3D>(trafoWorld),
                               worldVol, trVolArr, "World"));

    mtvpWorld->sign(GeometrySignature::Global);

    // Build and return tracking geometry
    return std::shared_ptr<TrackingGeometry>(
        new Acts::TrackingGeometry(mtvpWorld));
  }

  RotationMatrix3D rotation = RotationMatrix3D::Identity();
  std::shared_ptr<const RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const ISurfaceMaterial> surfaceMaterial = nullptr;

  std::vector<std::unique_ptr<const DetectorElementStub>> detectorStore = {};

  std::reference_wrapper<const GeometryContext> geoContext;
};
}  // namespace Test
}  // namespace Acts
