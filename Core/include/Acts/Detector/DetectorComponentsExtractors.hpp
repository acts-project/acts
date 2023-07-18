// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;
struct NavigationState;

/// Helper extractors: all portals
struct AllPortalsExtractor {
  /// Extract the portals from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  ///
  /// @return a vector of raw Portal pointers
  static const std::vector<const Portal*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState);
};

/// Helper extractors: all surfaces
struct AllSurfacesExtractor {
  /// Extract the surfaces from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices is an ignored index vector
  ///
  /// @return a vector of raw Surface pointers
  static const std::vector<const Surface*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      [[maybe_unused]] const std::vector<std::size_t>& indices = {});
};

/// Helper extractors: indexed surfaces
struct IndexedSurfacesExtractor {
  /// Extract the surfaces from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the surfaces store
  ///
  /// @note no out of boudns checking is done
  ///
  /// @return a vector of raw Surface pointers
  static const std::vector<const Surface*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState, const std::vector<size_t>& indices);
};

/// Helper extractors: all sub volumes of a volume
struct AllSubVolumesExtractor {
  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the volume store (ignored)
  ///
  /// @return a vector of raw DetectorVolume pointers
  static const std::vector<const DetectorVolume*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      [[maybe_unused]] const std::vector<std::size_t>& indices = {});
};

/// Helper extractors: indexed sub volume of a volume
struct IndexedSubVolumesExtractor {
  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the volume store
  ///
  /// @return a vector of raw DetectorVolume pointers
  static const std::vector<const DetectorVolume*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState, const std::vector<std::size_t>& indices);
};

}  // namespace Experimental
}  // namespace Acts
