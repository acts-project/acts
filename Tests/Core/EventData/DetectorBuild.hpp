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
#include "Acts/Detector/detail/DefaultDetectorElementBase.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
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

  ///
  /// @brief Stub implementation of the detector element
  ///
  class DetElem : public DetectorElementBase
  {
  public:
    /// @brief Constructor
    ///
    /// @param [in] trafo Transformation of the detector element
    /// @param [in] rBounds Rectangle boundaries of the plane surface
    /// @param [in] thickness Thickness of the detector element
    DetElem(std::shared_ptr<const Transform3D>     trafo,
            std::shared_ptr<const RectangleBounds> rBounds,
            double                                 thickness)
      : DetectorElementBase()
      , m_trafo(trafo)
      , m_surface(new PlaneSurface(rBounds, *this))
      , m_thickness(thickness)
    {
    }

    /// @brief Getter of the transformation
    virtual const Transform3D&
    transform() const
    {
      return *m_trafo;
    }

    /// @brief Getter of the surface
    virtual const Surface&
    surface() const
    {
      return *m_surface;
    }

    /// @brief Getter of the thickness
    virtual double
    thickness() const
    {
      return m_thickness;
    }

    // Pointer to the transformation
    std::shared_ptr<const Transform3D> m_trafo;
    // Surface related to the detector element
    Surface const* m_surface;
    // Thickness of the detector element
    double m_thickness;
  };

  /// @brief Builds a simple 4-layer detector with 2 pixel-like and 2
  /// double-strip-like detectors
  ///
  /// @return Pointer to the tracking geometry
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

    MaterialProperties matProp(
        352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
    std::shared_ptr<const SurfaceMaterial> surMat(
        new HomogeneousSurfaceMaterial(matProp));

    // Construct surfaces
    std::array<PlaneSurface*, 6> surfaces;
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];

      surfaces[i] = new PlaneSurface(
          rBounds,
          *(new DetElem(std::make_shared<const Transform3D>(trafo),
                        rBounds,
                        1. * units::_um)));
      surfaces[i]->setAssociatedMaterial(surMat);
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
      surfaces[i]->associateLayer(*layers[i]);
    }

    // Build volume for surfaces with negative x-values
    Transform3D trafoVol1(Transform3D::Identity());
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
    Transform3D trafoVol2(Transform3D::Identity());
    trafoVol2.translation() = Vector3D(1.5 * units::_m, 0., 0.);

    layVec.clear();
    for (i = 2; i < 6; i++) layVec.push_back(layers[i]);
    std::unique_ptr<const LayerArray> layArr2(
        layArrCreator.layerArray(layVec,
                                 1. * units::_m - 2. * units::_mm,
                                 2. * units::_m + 2. * units::_mm,
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
    trackVolume2->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                     trackVolume1,
                                     BoundarySurfaceFace::positiveFaceYZ);

    trackVolume1->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                     trackVolume2,
                                     BoundarySurfaceFace::negativeFaceYZ);

    // Build world volume
    Transform3D trafoWorld(Transform3D::Identity());
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
