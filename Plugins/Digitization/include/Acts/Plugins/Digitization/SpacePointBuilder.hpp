// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"

namespace Acts {

/// @struct SpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits on some detector elements need further treatment. This struct
/// serves as default structure of the process to take the digitized clusters on
/// a
/// detector element and provide the corresponding space point. The empty
/// class is used to forbid the usage of an arbitrary data type as template
/// parameter and enforces the implementation of explicit structures.
///
/// @note The choice of which kind of data should be treated in which way is
/// steered by the choice of the template parameter. This parameter represents a
/// structure that needs to store at least a cluster/multiple clusters and the
/// corresponding space point.
///
template <class S>
class SpacePointBuilder
{
};
}  // namespace Acts
