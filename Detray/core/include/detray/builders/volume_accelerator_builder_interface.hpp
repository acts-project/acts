// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <string_view>
#include <vector>

namespace detray {

/// @brief Interface for builders of volume acceleration structures.
///
/// A volume acceleration structure holds the daughter volumes of a "root"
/// volume. The root links to the acceleration structure through its
/// @c e_volume accelerator link. A root of @c dindex_invalid stands for the
/// top-level volumes (world volume lookup).
///
/// The structures can only be built after all volumes are built. The
/// @c detector_builder holds a collection of these builders and calls them in
/// the order in which they were added, after it built the volumes. Roots that
/// are not registered with any builder get a volume brute force accelerator.
///
/// The volume acceleration builders are independent of the volume builders
/// and of each other: Add only the builders that the detector needs.
///
/// @tparam detector_t the type of the detector under construction
/// @note The template parameter is not constrained, because the detector
/// declares some builders as friends while the detector type is incomplete.
template <typename detector_t>
class volume_accelerator_builder_interface {
 public:
  /// Default destructor
  virtual ~volume_accelerator_builder_interface() = default;

  /// @returns a name of the builder (used in log output)
  DETRAY_HOST
  virtual auto name() const -> std::string_view = 0;

  /// Register a root volume and its daughter volumes with the builder
  ///
  /// @param root index of the root volume (@c dindex_invalid for the
  ///             top-level volumes)
  /// @param daughters indices of the daughter volumes of the root
  DETRAY_HOST
  virtual void add_root(dindex root, const std::vector<dindex> &daughters) = 0;

  /// @returns the indices of the root volumes registered with the builder
  DETRAY_HOST
  virtual auto roots() const -> std::vector<dindex> = 0;

  /// Build the acceleration structures into the detector and set the
  /// accelerator links of the root volumes
  ///
  /// @param det the detector, which already contains all volumes
  DETRAY_HOST
  virtual void build(detector_t &det) = 0;
};

}  // namespace detray
