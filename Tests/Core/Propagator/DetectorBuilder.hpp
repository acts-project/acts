// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"

namespace Acts {
namespace Test {

  /// @brief Builds tracking geometry without anything in it
  ///
  /// @return Pointer to the tracking geometry
  std::shared_ptr<TrackingGeometry>
  buildVacDetector()
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

    // Build the transformation
    Transform3D trafoLay(Transform3D::Identity() * rotation);
    trafoLay.translation() = Vector3D(1. * units::_m, 0., 0.);

    // Build a dummy layer
    std::shared_ptr<const PlanarBounds> rBounds(
        new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    LayerPtr dummyLayer = PlaneLayer::create(
        std::make_shared<const Transform3D>(trafoLay), rBounds);

    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    std::unique_ptr<const LayerArray> layArr(
        layArrCreator.layerArray({dummyLayer},
                                 0.,
                                 2. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    // Build the volume
    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
        1. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    Transform3D trafoVac(Transform3D::Identity());
    trafoVac.translation() = Vector3D(1. * units::_m, 0., 0.);
    auto trackingVac
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVac),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArr),
                                 {},
                                 {},
                                 {},
                                 "Vacuum");

    return std::shared_ptr<TrackingGeometry>(new TrackingGeometry(trackingVac));
  }

  /// @brief Builds tracking geometry that contains one volume with material
  ///
  /// @return Pointer to the tracking geometry
  std::shared_ptr<TrackingGeometry>
  buildMatDetector()
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

    // Build the transformation
    Transform3D trafoLay(Transform3D::Identity() * rotation);
    trafoLay.translation() = Vector3D(1. * units::_m, 0., 0.);

    // Build a dummy layer
    std::shared_ptr<const PlanarBounds> rBounds(
        new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    LayerPtr dummyLayer = PlaneLayer::create(
        std::make_shared<const Transform3D>(trafoLay), rBounds);

    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    std::unique_ptr<const LayerArray> layArr(
        layArrCreator.layerArray({dummyLayer},
                                 0.,
                                 2. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    // Build the volume
    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
        1. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    std::shared_ptr<const Material> mat(
        new Material(352.8, 407., 9.012, 4., 1.848e-3));

    Transform3D trafoMat(Transform3D::Identity());
    trafoMat.translation() = Vector3D(1. * units::_m, 0., 0.);
    auto trackingMat
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoMat),
                                 boundsVol,
                                 mat,
                                 std::move(layArr),
                                 {},
                                 {},
                                 {},
                                 "Material");

    return std::shared_ptr<TrackingGeometry>(new TrackingGeometry(trackingMat));
  }

  std::shared_ptr<TrackingGeometry>
  buildVacMatVacDetector()
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

    // Build the transformation
    Transform3D trafoLay1(Transform3D::Identity() * rotation);
    trafoLay1.translation() = Vector3D(1. * units::_m, 0., 0.);
    Transform3D trafoLay2(Transform3D::Identity() * rotation);
    trafoLay2.translation() = Vector3D(3. * units::_m, 0., 0.);
    Transform3D trafoLay3(Transform3D::Identity() * rotation);
    trafoLay3.translation() = Vector3D(5. * units::_m, 0., 0.);

    // Build a dummy layer
    std::shared_ptr<const PlanarBounds> rBounds(
        new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    LayerPtr dummyLayer1 = PlaneLayer::create(
        std::make_shared<const Transform3D>(trafoLay1), rBounds);
    LayerPtr dummyLayer2 = PlaneLayer::create(
        std::make_shared<const Transform3D>(trafoLay2), rBounds);
    LayerPtr dummyLayer3 = PlaneLayer::create(
        std::make_shared<const Transform3D>(trafoLay3), rBounds);

    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    std::unique_ptr<const LayerArray> layArr1(
        layArrCreator.layerArray({dummyLayer1},
                                 0.,
                                 2. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binX));
    std::unique_ptr<const LayerArray> layArr2(
        layArrCreator.layerArray({dummyLayer2},
                                 2.,
                                 4. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binX));
    std::unique_ptr<const LayerArray> layArr3(
        layArrCreator.layerArray({dummyLayer3},
                                 4.,
                                 6. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    // Build the volume
    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
        1. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    Transform3D trafoVol1(Transform3D::Identity());
    trafoVol1.translation() = Vector3D(1. * units::_m, 0., 0.);
    Transform3D trafoVol2(Transform3D::Identity());
    trafoVol2.translation() = Vector3D(3. * units::_m, 0., 0.);
    Transform3D trafoVol3(Transform3D::Identity());
    trafoVol3.translation() = Vector3D(5. * units::_m, 0., 0.);

    std::shared_ptr<const Material> mat(
        new Material(352.8, 407., 9.012, 4., 1.848e-3));

    auto trackingVac1
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol1),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArr1),
                                 {},
                                 {},
                                 {},
                                 "Vacuum1");
    trackingVac1->sign(GeometrySignature::Global);
    auto trackingMat
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol2),
                                 boundsVol,
                                 mat,
                                 std::move(layArr2),
                                 {},
                                 {},
                                 {},
                                 "Material");
    trackingMat->sign(GeometrySignature::Global);
    auto trackingVac2
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol3),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArr3),
                                 {},
                                 {},
                                 {},
                                 "Vacuum2");
    trackingVac2->sign(GeometrySignature::Global);

    // Glue volumes
    trackingVac1->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                     trackingMat,
                                     BoundarySurfaceFace::negativeFaceYZ);
    trackingMat->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                    trackingVac1,
                                    BoundarySurfaceFace::positiveFaceYZ);
    trackingMat->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                    trackingVac2,
                                    BoundarySurfaceFace::negativeFaceYZ);
    trackingVac2->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                     trackingMat,
                                     BoundarySurfaceFace::positiveFaceYZ);

    // Build world volume
    Transform3D trafoWorld(Transform3D::Identity());
    trafoWorld.translation() = Vector3D(3. * units::_m, 0., 0.);

    auto worldVol = std::make_shared<const CuboidVolumeBounds>(
        3. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

    std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

    tapVec.push_back(
        std::make_pair(trackingVac1, Vector3D(1. * units::_m, 0., 0.)));
    tapVec.push_back(
        std::make_pair(trackingMat, Vector3D(3. * units::_m, 0., 0.)));
    tapVec.push_back(
        std::make_pair(trackingVac2, Vector3D(5. * units::_m, 0., 0.)));

    std::vector<double> binBoundaries
        = {0., 2. * units::_m, 4. * units::_m, 6. * units::_m};

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
