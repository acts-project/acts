// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryStatics.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H
#define ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H 1

#include "ACTS/Utilities/Definitions.hpp"

/// Define statics for Geometry in Tracking
/// 
namespace Acts {

// transformations

static Transform3D s_idTransform
    = Transform3D::Identity();  //!< idendity transformation
static Rotation3D s_idRotation
    = Acts::Rotation3D::Identity();  //!< idendity rotation

// axis system
static Vector3D s_xAxis(1, 0, 0);  //!< global x Axis;
static Vector3D s_yAxis(0, 1, 0);  //!< global y Axis;
static Vector3D s_zAxis(0, 0, 1);  //!< global z Axis;

// unit vectors
static Vector2D s_origin2D(0.,0.);

// origin

static Vector3D s_origin(0, 0, 0);  //!< origin position

static double helper[9] = {0., 1., 0., 1., 0., 0., 0., 0., -1.};

static Acts::RotationMatrix3D s_idRotationZinverse(helper);
}

#endif  // ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H
