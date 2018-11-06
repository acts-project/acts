// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoxGeometryBuilderTest

#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Tools/BoxGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

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

  BOOST_AUTO_TEST_CASE(BoxGeometryBuilderTest)
  {

    // Construct builder
    BoxGeometryBuilder bgb;

    // Create configurations for surfaces
    std::vector<BoxGeometryBuilder::SurfaceConfig> surfaceConfig;
    for (unsigned int i = 1; i < 5; i++) {
      // Position of the surfaces
      BoxGeometryBuilder::SurfaceConfig cfg;
      cfg.position = {i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;

      surfaceConfig.push_back(cfg);
    }

    // Test that there are actually 4 surface configurations
    BOOST_TEST(surfaceConfig.size() == 4);

    // Test that 4 sensitive surfaces can be built
    for (const auto& cfg : surfaceConfig) {
      PlaneSurface* pSur = bgb.buildSurface<DetElem>(cfg);
      BOOST_TEST(pSur != nullptr);
      BOOST_TEST(pSur->center() == cfg.position);
      BOOST_TEST(pSur->associatedMaterial() != nullptr);
      BOOST_TEST(pSur->associatedDetectorElement() != nullptr);
    }

    // Test that 4 passive surfaces can be built
    for (const auto& cfg : surfaceConfig) {
      PlaneSurface* pSur = bgb.buildSurface<>(cfg);
      BOOST_TEST(pSur != nullptr);
      BOOST_TEST(pSur->center() == cfg.position);
      BOOST_TEST(pSur->associatedMaterial() != nullptr);
      BOOST_TEST(pSur->associatedDetectorElement() == nullptr);
    }

    ////////////////////////////////////////////////////////////////////
    // Build layer configurations
    std::vector<BoxGeometryBuilder::LayerConfig> layerConfig;
    for (auto& sCfg : surfaceConfig) {
      BoxGeometryBuilder::LayerConfig cfg;
      cfg.layerThickness = 1. * units::_mm;
      cfg.surfaceCfg     = sCfg;
      layerConfig.push_back(cfg);
    }

    // Test that there are actually 4 layer configurations
    BOOST_TEST(layerConfig.size() == 4);

    // Test that 4 layers with sensitive surfaces can be built
    for (auto& cfg : layerConfig) {
      LayerPtr layer = bgb.buildLayer<DetElem>(cfg);
      BOOST_TEST(layer != nullptr);
      BOOST_TEST(layer->thickness() == 1. * units::_mm);
      BOOST_TEST(cfg.surface != nullptr);
      BOOST_TEST(layer->surfaceArray()->size() == 1);
      BOOST_TEST(layer->layerType() == LayerType::active);
    }

    // Test that 4 layers with passive surfaces can be built
    for (auto& cfg : layerConfig) {
      LayerPtr layer = bgb.buildLayer<>(cfg);
      BOOST_TEST(layer != nullptr);
      BOOST_TEST(layer->thickness() == 1. * units::_mm);
      BOOST_TEST(cfg.surface != nullptr);
      BOOST_TEST(layer->surfaceArray()->size() == 1);
      BOOST_TEST(layer->layerType() == LayerType::passive);
    }

    ////////////////////////////////////////////////////////////////////
    // Build volume configuration
    BoxGeometryBuilder::VolumeConfig volumeConfig;
    volumeConfig.position = {2.5 * units::_m, 0., 0.};
    volumeConfig.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig.layerCfg = layerConfig;
    volumeConfig.binningValue = BinningValue::binX;
    volumeConfig.name         = "Test volume";

    // Test the building
    std::shared_ptr<TrackingVolume> trVol = bgb.buildVolume(volumeConfig);
    BOOST_TEST(volumeConfig.layers.size() == 4);
    BOOST_TEST(trVol->confinedLayers()->arrayObjects().size()
               == volumeConfig.layers.size() * 2
                   + 1);  // #layers = navigation + material layers
    BOOST_TEST(trVol->volumeName() == volumeConfig.name);

    ////////////////////////////////////////////////////////////////////
    // Build TrackingGeometry configuration

    std::vector<BoxGeometryBuilder::SurfaceConfig> surfaceConfig2;
    for (unsigned int i = 1; i < 5; i++) {
      // Position of the surfaces
      BoxGeometryBuilder::SurfaceConfig cfg;
      cfg.position = {-i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;

      surfaceConfig2.push_back(cfg);
    }

    std::vector<BoxGeometryBuilder::LayerConfig> layerConfig2;
    for (auto& sCfg : surfaceConfig2) {
      BoxGeometryBuilder::LayerConfig cfg;
      cfg.layerThickness = 1. * units::_mm;
      cfg.surfaceCfg     = sCfg;
      layerConfig2.push_back(cfg);
    }
    BoxGeometryBuilder::VolumeConfig volumeConfig2;
    volumeConfig2.position = {-2.5 * units::_m, 0., 0.};
    volumeConfig2.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig2.layerCfg = layerConfig;
    volumeConfig2.binningValue = BinningValue::binX;
    volumeConfig2.name         = "Test volume2";

    BoxGeometryBuilder::Config config;
    config.position  = {0., 0., 0.};
    config.length    = {10. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {volumeConfig, volumeConfig2};

    bgb.buildTrackingGeometry(config);
  }
}  // namespace Test
}  // namespace Acts
