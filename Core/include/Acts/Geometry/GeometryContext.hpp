// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Set the Geometry Context PLUGIN
#ifdef ACTS_CORE_GEOMETRYCONTEXT_PLUGIN
#include ACTS_CORE_GEOMETRYCONTEXT_PLUGIN
#else

#include "Acts/Utilities/detail/ContextType.hpp"

namespace Acts {

/// @brief This is the central definition of the Acts
/// payload object regarding detector geometry status (e.g. alignment)
///
/// It is propagated through the code to allow for event/thread
/// dependent geometry changes

using GeometryContext = ContextType;

}  // namespace Acts

#endif
