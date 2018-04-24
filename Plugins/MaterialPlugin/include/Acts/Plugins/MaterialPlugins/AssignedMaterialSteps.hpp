// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AssignedMaterialSteps.h, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"

namespace Acts {

/// This struct is used to assign a AssignedMaterialSteps to a mapped position
///
struct AssignedMaterialSteps
{

  GeometryID assignedGeoID;     ///< this is the geo ID of the assigned surface
  Vector3D   assignedPosition;  ///< this is the position of intersection
  std::vector<MaterialStep> assignedSteps;  ///< this is step information

  // simple constructor
  AssignedMaterialSteps(GeometryID      geoID    = GeometryID(),
                        const Vector3D& position = Vector3D(0., 0., 0),
                        const std::vector<MaterialStep>& steps = {})
    : assignedGeoID(geoID), assignedPosition(position), assignedSteps(steps)
  {
  }
};
}