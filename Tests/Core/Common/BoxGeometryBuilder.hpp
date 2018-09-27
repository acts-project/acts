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

	class BoxGeometryBuilder
	{
		public:
			BoxGeometryBuilder();
			
			template<typename DetectorElement_t>
			std::shared_ptr<TrackingGeometry>
			buildGeometry() const;
			
			template<typename DetectorElement_t>
			std::shared_ptr<TrackingGeometry>
			buildGeometry(const std::vector<Vector3D>& pixelSurfaces, const std::vector<Vector3D>& stripSurfaces, const std::pair<double, double>& detectorLength) const;
			
			template<typename DetectorElement_t>
			std::shared_ptr<TrackingGeometry>
			buildGeometry(const std::vector<Vector3D>& pixelSurfaces, const std::vector<Vector3D>& stripSurfaces, const double distStrips, const std::pair<double, double>& detectorLength) const;
		
		private:
		
			template<typename DetectorElement_t>
			std::vector<PlaneSurface*>
			buildSurfaces(const std::vector<Vector3D>& surfacePos) const;
			
			RotationMatrix3D rotation;
			std::shared_ptr<const RectangleBounds> rBounds;
			std::shared_ptr<const SurfaceMaterial> surMat;
	};

	BoxGeometryBuilder::BoxGeometryBuilder()
	{
		// Construct the rotation
		double           rotationAngle = M_PI * 0.5;
		Vector3D         xPos(cos(rotationAngle), 0., sin(rotationAngle));
		Vector3D         yPos(0., 1., 0.);
		Vector3D         zPos(-sin(rotationAngle), 0., cos(rotationAngle));
		rotation.col(0) = xPos;
		rotation.col(1) = yPos;
		rotation.col(2) = zPos;
		
		// Boundaries of the surfaces
		rBounds = std::make_shared<const RectangleBounds>(
			RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

		// Material of the surfaces
		MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
		surMat = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(matProp));
	}

  /// @brief Builds a simple 4-layer detector with 2 pixel-like and 2
  /// double-strip-like detectors
  ///
  /// @return Pointer to the tracking geometry
  template<typename DetectorElement_t>
  std::shared_ptr<TrackingGeometry>
  BoxGeometryBuilder::buildGeometry() const
  {
    // Set translation vectors
    double                eps = 1. * units::_mm;
    std::vector<Vector3D> translations;
    translations.push_back({-2. * units::_m, 0., 0.});
    translations.push_back({-1. * units::_m, 0., 0.});
    translations.push_back({1. * units::_m - eps, 0., 0.});
    translations.push_back({1. * units::_m + eps, 0., 0.});
    translations.push_back({2. * units::_m - eps, 0., 0.});
    translations.push_back({2. * units::_m + eps, 0., 0.});

    // Construct surfaces
    std::array<PlaneSurface*, 6> surfaces;
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = translations[i];

      surfaces[i] = new PlaneSurface(
          rBounds,
          *(new DetectorElement_t(std::make_shared<const Transform3D>(trafo),
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

  /// @brief Builds a simple 4-layer detector with 2 pixel-like and 2
  /// double-strip-like detectors
  ///
  /// @return Pointer to the tracking geometry
  template<typename DetectorElement_t>
  std::shared_ptr<TrackingGeometry>
  BoxGeometryBuilder::buildGeometry(const std::vector<Vector3D>& pixelSurfaces, const std::vector<Vector3D>& stripSurfaces, const std::pair<double, double>& detectorLength) const
  {
    // Construct surfaces
    std::vector<PlaneSurface*> surfacesP = buildSurfaces<DetectorElement_t>(pixelSurfaces);
    std::vector<PlaneSurface*> surfacesS = buildSurfaces<DetectorElement_t>(stripSurfaces);

    // Construct layers
    std::vector<LayerPtr> layersP, layersS;
    layersP.reserve(pixelSurfaces.size());
    layersS.reserve(stripSurfaces.size());
    for (unsigned int i = 0; i < pixelSurfaces.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = pixelSurfaces[i];

      std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfacesP[i]));

      layersP.push_back(PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
                                     rBounds,
                                     std::move(surArray),
                                     1. * units::_mm));
      surfacesP[i]->associateLayer(*layersP[i]);
    }
    for (unsigned int i = 0; i < stripSurfaces.size(); i++) {
      Transform3D trafo(Transform3D::Identity() * rotation);
      trafo.translation() = stripSurfaces[i];

      std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfacesS[i]));

      layersS.push_back(PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
                                     rBounds,
                                     std::move(surArray),
                                     1. * units::_mm));
      surfacesS[i]->associateLayer(*layersS[i]);
    }

    // Build volume for surfaces with negative x-values
    Transform3D trafoVolP(Transform3D::Identity());
    trafoVolP.translation() = Vector3D(detectorLength.first * 0.5, 0., 0.);

    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
        detectorLength.first, 0.5 * units::_m, 0.5 * units::_m);

    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    LayerVector layVec;
    for(unsigned int i = 0; i < layersP.size(); i++)
		layVec.push_back(layersP[i]);

    std::unique_ptr<const LayerArray> layArrP(
        layArrCreator.layerArray(layVec,
                                 pixelSurfaces[0].x() - 1. * units::_mm,
                                 pixelSurfaces.back().x() + 1. * units::_mm,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    auto trackVolumeP
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVolP),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArrP),
                                 layVec,
                                 {},
                                 {},
                                 "Volume 1");
    trackVolumeP->sign(GeometrySignature::Global);

    // Build volume for surfaces with positive x-values
    Transform3D trafoVolS(Transform3D::Identity());
    trafoVolS.translation() = Vector3D(detectorLength.second * 0.5, 0., 0.);

    layVec.clear();
    for(unsigned int i = 0; i < layersS.size(); i++)
		layVec.push_back(layersS[i]);
		
    std::unique_ptr<const LayerArray> layArrS(
        layArrCreator.layerArray(layVec,
                                 stripSurfaces[0].x() - 1. * units::_mm,
                                 stripSurfaces.back().x() + 1. * units::_mm,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    auto trackVolumeS
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVolS),
                                 boundsVol,
                                 nullptr,
                                 std::move(layArrS),
                                 layVec,
                                 {},
                                 {},
                                 "Volume 2");
    trackVolumeS->sign(GeometrySignature::Global);

    // Glue volumes
    trackVolumeS->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                     trackVolumeP,
                                     BoundarySurfaceFace::positiveFaceYZ);

    trackVolumeP->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                     trackVolumeS,
                                     BoundarySurfaceFace::negativeFaceYZ);

    // Build world volume
    Transform3D trafoWorld(Transform3D::Identity());
    trafoWorld.translation() = Vector3D(0., 0., 0.);

    auto worldVol = std::make_shared<const CuboidVolumeBounds>(
        (detectorLength.first + detectorLength.second) * 0.5, 0.5 * units::_m, 0.5 * units::_m);

    std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

    tapVec.push_back(
        std::make_pair(trackVolumeP, Vector3D(detectorLength.first * 0.5, 0., 0.)));
    tapVec.push_back(
        std::make_pair(trackVolumeS, Vector3D(detectorLength.second * 0.5, 0., 0.)));

    std::vector<double> binBoundaries = {detectorLength.first, 0., detectorLength.second};

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

	template<typename DetectorElement_t>
	std::vector<PlaneSurface*>
	BoxGeometryBuilder::buildSurfaces(const std::vector<Vector3D>& surfacePos) const
	{
		// Construct surfaces
		std::vector<PlaneSurface*> surfaces;
		surfaces.reserve(surfacePos.size());
		
		for (unsigned int i = 0; i < surfacePos.size(); i++) {
		  Transform3D trafo(Transform3D::Identity() * rotation);
		  trafo.translation() = surfacePos[i];

		  surfaces.push_back(new PlaneSurface(
			  rBounds,
			  *(new DetectorElement_t(std::make_shared<const Transform3D>(trafo),
							rBounds,
							1. * units::_um))));
		  surfaces[i]->setAssociatedMaterial(surMat);
		}
		return surfaces;
	}
}  // namespace Test
}  // namespace Acts
