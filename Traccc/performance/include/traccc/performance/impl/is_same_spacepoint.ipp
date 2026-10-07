// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/performance/details/is_same_scalar.hpp"

// Project include(s).
#include "traccc/edm/spacepoint_collection.hpp"

namespace traccc::details {

/// @c traccc::is_same_object specialisation for @c traccc::edm::spacepoint
template <typename T>
class is_same_object<edm::spacepoint<T>> {
 public:
  /// Constructor with a reference object, and an allowed uncertainty
  is_same_object(const edm::spacepoint<T>& ref, scalar unc = float_epsilon)
      : m_ref(ref), m_unc(unc) {}

  /// Specialised implementation for @c traccc::spacepoint
  bool operator()(const edm::spacepoint<T>& obj) const {
    return (is_same_scalar(obj.x(), m_ref.x(), m_unc) &&
            is_same_scalar(obj.y(), m_ref.y(), m_unc) &&
            is_same_scalar(obj.z(), m_ref.z(), m_unc));
  }

 private:
  /// The reference object
  const edm::spacepoint<T> m_ref;
  /// The uncertainty
  scalar m_unc;

};  // class is_same_object<spacepoint>

}  // namespace traccc::details
