// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"

// Local include(s).
#include "traccc/performance/details/is_same_object.hpp"

namespace traccc::details {

/// @c traccc::details::comparator_factory specialisation for
/// @c traccc::edm::seed
template <typename T>
class comparator_factory<edm::seed<T>> {
 public:
  /// Constructor with all necessary arguments
  comparator_factory(
      const edm::spacepoint_collection::const_view& ref_spacepoints,
      const edm::spacepoint_collection::const_view& test_spacepoints)
      : m_ref_spacepoints(ref_spacepoints),
        m_test_spacepoints(test_spacepoints) {}

  /// Instantiate an instance of a comparator object
  is_same_object<edm::seed<T>> make_comparator(
      const edm::seed<T>& ref, scalar unc = float_epsilon) const {
    return is_same_object<edm::seed<T>>(m_ref_spacepoints, m_test_spacepoints,
                                        ref, unc);
  }

 private:
  /// Spacepoint container for the reference seeds
  const edm::spacepoint_collection::const_view m_ref_spacepoints;
  /// Spacepoint container for the test seeds
  const edm::spacepoint_collection::const_view m_test_spacepoints;

};  // class comparator_factory

}  // namespace traccc::details
