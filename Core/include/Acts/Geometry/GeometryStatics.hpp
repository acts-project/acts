// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"

/// Define statics for Geometry in Tracking
///
namespace Acts {

// Transformations

static const Transform3D s_idTransform =
    Transform3D::Identity();  //!< idendity transformation

}  // namespace Acts
