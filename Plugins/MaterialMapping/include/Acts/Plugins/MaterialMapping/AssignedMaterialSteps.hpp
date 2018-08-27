// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AssignedMaterialSteps.h, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Plugins/MaterialMapping/MaterialStep.hpp"
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
  AssignedMaterialSteps(GeometryID geoID                = GeometryID(),
                        Vector3D   position             = Vector3D(0., 0., 0),
                        std::vector<MaterialStep> steps = {})
    : assignedGeoID(std::move(geoID))
    , assignedPosition(std::move(position))
    , assignedSteps(std::move(steps))
  {
  }
};
}