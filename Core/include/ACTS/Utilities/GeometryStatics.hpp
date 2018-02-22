// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
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

static const Transform3D s_idTransform
    = Transform3D::Identity();  //!< idendity transformation
static const Rotation3D s_idRotation
    = Acts::Rotation3D::Identity();  //!< idendity rotation

// axis system
static const Vector3D s_xAxis(1, 0, 0);  //!< global x Axis;
static const Vector3D s_yAxis(0, 1, 0);  //!< global y Axis;
static const Vector3D s_zAxis(0, 0, 1);  //!< global z Axis;

// unit vectors
static const Vector2D s_origin2D(0., 0.);

// origin

static const Vector3D s_origin(0, 0, 0);  //!< origin position

static const double helper[9] = {0., 1., 0., 1., 0., 0., 0., 0., -1.};

static const Acts::RotationMatrix3D s_idRotationZinverse(helper);
}

#endif  // ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H
