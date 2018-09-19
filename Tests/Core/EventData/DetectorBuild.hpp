// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"

namespace Acts {
namespace Test {

  std::shared_ptr<TrackingGeometry>
  buildGeometry()
  {
    // Construct the rotation
    RotationMatrix3D rotation;
    double           rotationAngle = M_PI * 0.5;
    Vector3D         xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D         yPos(0., 1., 0.);
    Vector3D         zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Set translation vectors
    double                eps = 1. * units::_mm;
    std::vector<Vector3D> translations;
    translations.push_back({-2. * units::_m, 0., 0.});
    translations.push_back({-1. * units::_m, 0., 0.});
    translations.push_back({1. * units::_m - eps, 0., 0.});
    translations.push_back({1. * units::_m + eps, 0., 0.});
    translations.push_back({2. * units::_m - eps, 0., 0.});
    translations.push_back({2. * units::_m + eps, 0., 0.});

    // Boundaries of the surface
    std::shared_ptr<const RectangleBounds> rBounds(
        new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
	
    // Construct surfaces
    std::array<PlaneSurface*, 6> surfaces;
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];
      surfaces[i] = new PlaneSurface(std::make_shared<const Transform3D>(trafo),
                                     rBounds);
    }

    // Construct layers
    std::array<LayerPtr, 6> layers;
    for (i = 0; i < 6; i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];

      std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfaces[i]));

      layers[i] = PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
                                     rBounds,
                                     std::move(surArray),
                                     1. * units::_mm);
    }

    // Build volume for surfaces with negative x-values
    Transform3D trafoVol1(
        Transform3D::Identity()
        * rotation);  // TODO: is rotation and the cube volume a duplication?
    trafoVol1.translation() = Vector3D(-1.5 * units::_m, 0., 0.);

    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
        1.5 * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    LayerVector layVec;
    layVec.push_back(layers[0]);
    layVec.push_back(layers[1]);
    std::unique_ptr<const LayerArray> layArr1(
        layArrCreator.layerArray(layVec,
                                 -2. * units::_m - 1. * units::_mm,
                                 -1. * units::_m + 1. * units::_mm,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    auto trackVolume1
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol1),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArr1),
                                 layVec,
                                 {},
                                 {},
                                 "Volume 1");
    trackVolume1->sign(GeometrySignature::Global);	
	
    // Build volume for surfaces with positive x-values
    Transform3D trafoVol2(Transform3D::Identity() * rotation);
    trafoVol2.translation() = Vector3D(1.5 * units::_m, 0., 0.);

    layVec.clear();
    for (i = 2; i < 6; i++) layVec.push_back(layers[i]);
    std::unique_ptr<const LayerArray> layArr2(
        layArrCreator.layerArray(layVec,
                                 1. * units::_m - 1. * units::_mm,
                                 2. * units::_m + 1. * units::_mm,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    auto trackVolume2
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol2),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArr2),
                                 layVec,
                                 {},
                                 {},
                                 "Volume 2");
    trackVolume2->sign(GeometrySignature::Global);

    // Glue volumes
    trackVolume1->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                     trackVolume2,
                                     BoundarySurfaceFace::negativeFaceXY);

    trackVolume2->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                     trackVolume1,
                                     BoundarySurfaceFace::positiveFaceXY);

    // Build world volume
    Transform3D trafoWorld;
    trafoWorld.translation() = Vector3D(0., 0., 0.);

    auto worldVol = std::make_shared<const CuboidVolumeBounds>(
        3. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

    tapVec.push_back(
        std::make_pair(trackVolume1, Vector3D(-1.5 * units::_m, 0., 0.)));
    tapVec.push_back(
        std::make_pair(trackVolume2, Vector3D(1.5 * units::_m, 0., 0.)));

    std::vector<double> binBoundaries = {-3. * units::_m, 0., 3. * units::_m};

    BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
    std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

    std::shared_ptr<const TrackingVolumeArray> trVolArr(
        new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

    MutableTrackingVolumePtr mtvpWorld(
        TrackingVolume::create(std::make_shared<const Transform3D>(trafoWorld),
                               worldVol,
                               trVolArr,
                               "World"));

    mtvpWorld->sign(GeometrySignature::Global);

    // Build and return tracking geometry
    return std::shared_ptr<TrackingGeometry>(
        new Acts::TrackingGeometry(mtvpWorld));
  }

}  // namespace Test
}  // namespace Acts
