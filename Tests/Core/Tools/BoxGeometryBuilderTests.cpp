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
#include "Acts/Tools/BoxGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Material/Material.hpp"

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
	std::vector<BoxGeometryBuilder::Config> config;
	for(unsigned int i = 1; i < 5; i++)
	{
		// Position of the surfaces
		BoxGeometryBuilder::Config cfg;
		cfg.surfaceCfg.position = {i * units::_m, 0., 0.};
				
		// Rotation of the surfaces
		double           rotationAngle = M_PI * 0.5;
		Vector3D         xPos(cos(rotationAngle), 0., sin(rotationAngle));
		Vector3D         yPos(0., 1., 0.);
		Vector3D         zPos(-sin(rotationAngle), 0., cos(rotationAngle));
		cfg.surfaceCfg.rotation.col(0) = xPos;
		cfg.surfaceCfg.rotation.col(1) = yPos;
		cfg.surfaceCfg.rotation.col(2) = zPos;
		
		// Boundaries of the surfaces
		cfg.surfaceCfg.rBounds = std::make_shared<const RectangleBounds>(
			RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

		// Material of the surfaces
		MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
		cfg.surfaceCfg.surMat = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(matProp));
		
		// Thickness of the detector element
		cfg.surfaceCfg.thickness = 1. * units::_um;
		
		config.push_back(cfg);
	}
	
   BOOST_TEST(config.size() == 4);
   
   bgb.buildSensitiveSurfaces<DetElem>(config);
   for(const auto& cfg : config)
	BOOST_TEST(cfg.surface != nullptr);
   bgb.buildPassiveSurfaces(config);
   for(const auto& cfg : config)
	BOOST_TEST(cfg.surface != nullptr);
	
	bgb.buildLayers(config);
	for(const auto& cfg : config)
		BOOST_TEST(cfg.layer != nullptr);
  }
}  // namespace Test
}  // namespace Acts
