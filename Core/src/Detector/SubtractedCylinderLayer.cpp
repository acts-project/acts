// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedCylinderLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/SubtractedCylinderLayer.hpp"

Acts::SubtractedCylinderLayer::SubtractedCylinderLayer(const Acts::SubtractedCylinderSurface* subCyl,
						                              double thickness,
						                              int laytyp) :
  SubtractedCylinderSurface(*subCyl),
  Layer(nullptr, thickness, nullptr, nullptr, laytyp)
{}


Acts::SubtractedCylinderLayer::SubtractedCylinderLayer(const Acts::SubtractedCylinderLayer& clay, const Acts::Transform3D& transf):
  SubtractedCylinderSurface(clay, transf),
  Layer(clay)
{}
    
const Acts::SubtractedCylinderSurface& Acts::SubtractedCylinderLayer::surfaceRepresentation() const
{
  return (*this);
}
