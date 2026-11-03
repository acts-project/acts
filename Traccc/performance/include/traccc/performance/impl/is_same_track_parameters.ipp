// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Library include(s).
#include "traccc/edm/track_parameters.hpp"

namespace traccc::details {

/// @c traccc::is_same_object specialisation for
/// @c traccc::bound_track_parameters
template <>
class is_same_object<bound_track_parameters<>> {
 public:
  /// Constructor with a reference object, and an allowed uncertainty
  is_same_object(const bound_track_parameters<>& ref,
                 scalar unc = float_epsilon);

  /// Specialised implementation for @c traccc::bound_track_parameters
  bool operator()(const bound_track_parameters<>& obj) const;

 private:
  /// The reference object
  std::reference_wrapper<const bound_track_parameters<>> m_ref;
  /// The uncertainty
  scalar m_unc;

};  // class is_same_object<bound_track_parameters>

}  // namespace traccc::details
