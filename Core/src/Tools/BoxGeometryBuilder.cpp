// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoxGeometryBuilder.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/BoxGeometryBuilder.hpp"

void
Acts::BoxGeometryBuilder::buildPassiveSurfaces(std::vector<Config>& config) const
{
	PlaneSurface* surface;
	
	for(Config& cfg : config)
	{
		// Build transformation
	  Transform3D trafo(Transform3D::Identity() * cfg.surfaceCfg.rotation);
	  trafo.translation() = cfg.surfaceCfg.position;

		// Create and store surface
	  surface = new PlaneSurface(std::make_shared<const Transform3D>(trafo), cfg.surfaceCfg.rBounds);
	  surface->setAssociatedMaterial(cfg.surfaceCfg.surMat);
	  cfg.surface = surface;
	}
}
			
void
Acts::BoxGeometryBuilder::buildLayers(std::vector<Config>& config) const
{
	// TODO: thickness of layer
	for (Config& cfg : config) {
		// Build transformation centered at the surface position
	  Transform3D trafo(Transform3D::Identity() * cfg.surfaceCfg.rotation);
	  trafo.translation() = cfg.surfaceCfg.position;
	  
	  // Get (1!) surface, build layer and store it
	  std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(cfg.surface));

	  cfg.layer = PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
									 cfg.surfaceCfg.rBounds,
									 std::move(surArray),
									 1. * units::_mm);
	  cfg.surface->associateLayer(*cfg.layer);
	}
}
