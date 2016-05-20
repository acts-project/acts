// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedPlaneLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/SubtractedPlaneLayer.hpp"

Acts::SubtractedPlaneLayer::SubtractedPlaneLayer(const SubtractedPlaneSurface* subtrPlaneSurf,
                                                double thickness,
                                                Acts::OverlapDescriptor* olap,
                                                int laytyp) :
  SubtractedPlaneSurface(*subtrPlaneSurf),
  Layer(nullptr, thickness, olap, nullptr, laytyp) 
{}
  
Acts::SubtractedPlaneLayer::SubtractedPlaneLayer(const Acts::SubtractedPlaneLayer& play, const Acts::Transform3D& transf):
  SubtractedPlaneSurface(play,transf),
  Layer(play)
{}

const Acts::SubtractedPlaneSurface& Acts::SubtractedPlaneLayer::surfaceRepresentation() const
{
  return (*this);
}
