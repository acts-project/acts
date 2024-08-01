// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Set the Mangetic Field Context PLUGIN
#ifdef ACTS_CORE_MAGFIELDCONTEXT_PLUGIN
#include ACTS_CORE_MAGFIELDCONTEXT_PLUGIN
#else

#include "Acts/Utilities/detail/ContextType.hpp"

namespace Acts {

/// @brief This is the central definition of the Acts
/// payload object regarding magnetic field status
///
/// It is propagated through the code to allow for event/thread
/// dependent magnetic field changes

using MagneticFieldContext = ContextType;

}  // namespace Acts

#endif
