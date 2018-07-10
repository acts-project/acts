// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

/// @class DiscBounds
///
/// common base class for all bounds that are in a r/phi frame
///  - simply introduced to avoid wrong bound assigments to surfaces

class DiscBounds : public SurfaceBounds
{
};

}  // namespace