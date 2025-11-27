// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Set the Mangetic Field Context PLUGIN
#ifdef ACTS_CORE_MAGFIELDCONTEXT_PLUGIN
#include ACTS_CORE_MAGFIELDCONTEXT_PLUGIN
#else

#include "Acts/Utilities/detail/ContextType.hpp"

namespace Acts {

/// @ingroup context magnetic_field
/// @brief Context object for lookup of magnetic field values
///
/// The magnetic field context is an opaque type which contains experiment
/// specific event context information. This can be used to supply event
/// dependent data to the magnetic field instance, in case it is needed to
/// provide correct field values. The library itself does not make any
/// assumptions on the content of this context type (it is implemented using
/// `std::any`), but passes a reference through the call-chain to the field
/// implementation. An experiment specific field implementation is then expected
/// to performa cast to the concrete type, and use the contents.
///
/// An example use case of the context could be to look up conditions data /
/// records for the value of the magnetic field at the time of the event.
class MagneticFieldContext : public ContextType {
 public:
  /// Inherit all constructors
  using ContextType::ContextType;
  using ContextType::operator=;
};

}  // namespace Acts

#endif
