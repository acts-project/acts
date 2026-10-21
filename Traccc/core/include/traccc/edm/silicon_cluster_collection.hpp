// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/edm/container.hpp>

namespace traccc::edm {

/// Interface for the @c traccc::edm::silicon_cluster_collection class.
///
/// It provides the API that users would interact with, while using the
/// column(s)/array(s) defined in @c traccc::edm::silicon_cluster_collection.
///
template <typename BASE>
class silicon_cluster : public BASE {
 public:
  /// Inherit the base class's constructor(s)
  using BASE::BASE;

  /// @name Cell Information
  /// @{

  /// Indices of the silicon cells belonging to a cluster (non-const)
  ///
  /// @return A (non-const) jagged vector of <tt>unsigned int</tt> values
  ///
  TRACCC_HOST_DEVICE
  auto& cell_indices() { return BASE::template get<0>(); }
  /// Indices of the silicon cells belonging to a cluster (const)
  ///
  /// @return A (const) jagged vector of <tt>unsigned int</tt> values
  ///
  TRACCC_HOST_DEVICE
  const auto& cell_indices() const { return BASE::template get<0>(); }

};  // class silicon_cell_collection_interface

/// SoA container describing reconstructed silicon clusters
using silicon_cluster_collection =
    vecmem::edm::container<silicon_cluster,
                           vecmem::edm::type::jagged_vector<unsigned int> >;

}  // namespace traccc::edm
