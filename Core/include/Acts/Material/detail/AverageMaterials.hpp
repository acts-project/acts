// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"

namespace Acts {
namespace detail {

/// Compute the average properties for a combined slab of two materials.
///
/// @params slab1 Properties of the first material slab
/// @params slab2 Properties of the second material slab
/// @returns Material slab with the combined thickness and average parameters
///
/// The averaged material slab has the combined thickness of the two input slabs
/// and assumes the two input materials are homogeneously and continuously mixed
/// throughout the slab.
MaterialSlab combineSlabs(const MaterialSlab& slab1, const MaterialSlab& slab2);

}  // namespace detail
}  // namespace Acts
