// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <memory>
#include <ostream>
#include <utility>
#include <vector>

namespace Acts {

class Surface;

/// The detector blueprint
class LayerBlueprint {
 public:
  /// Constructor with all/part of the surfaces
  ///
  /// @param gctx The geometry context
  /// @param restriction The extent restriction against which this is checked
  /// @param envelope The envelope in all bin directions
  /// @param surfaces The surfaces in to be contained by this layer volume
  LayerBlueprint(
      const GeometryContext& gctx, const Extent& restriction = Extent(true),
      const std::vector<std::pair<ActsScalar, ActsScalar>>& envelope =
          std::vector<std::pair<ActsScalar, ActsScalar>>(size_t(binValues),
                                                         {0., 0.}),
      const std::vector<std::shared_ptr<Surface>>& surfaces = {});

  LayerBlueprint() = default;

  /// Add a surface, this will also check the restrictive and increase
  /// actual the extent
  ///
  /// If it fails the restriction, it will not be added
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surface The surface which is added to the layer BluePrint
  void add(const GeometryContext& gctx, std::shared_ptr<Surface> surface);

  /// Get the parameters : min
  /// @param bval The accessed binning value
  /// @param addenv The steering if enevlope is added or not
  ActsScalar min(BinningValue bval, bool addenv = true) const;

  /// Get the  parameters : max
  /// @param bval The accessed binning value
  /// @param addenv The steering if enevlope is added or not
  ActsScalar max(BinningValue bval, bool addenv = true) const;

  /// Get the  parameters : max
  /// @param bval The accessed binning value
  /// @param addenv The steering if enevlope is added or not
  ActsScalar medium(BinningValue bval, bool addenv = true) const;

  /// Get the  parameters : max
  /// @param bval The accessed binning value
  /// @param addenv The steering if enevlope is added or not
  ActsScalar range(BinningValue bval, bool addenv = true) const;

  /// @return the current extent
  const Extent& extent() const;

  /// @return the vector of surfaces
  const std::vector<std::shared_ptr<Surface>>& surfaces() const;

  /// Output to ostream
  /// @param sl the input ostream
  std::ostream& toStream(std::ostream& sl) const;

 private:
  /// Helper method which performs the actual min/max calculation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The surfaces to build this protolayer out of
  void measure(const GeometryContext& gctx,
               const std::vector<std::shared_ptr<Surface>>& surfaces);

  /// The restriction extent (surfaces are only accepted if they fall within)
  Extent m_restriction;

  /// The actual extent (measured from surfaces)
  Extent m_extent;

  /// Envelope
  std::vector<std::pair<ActsScalar, ActsScalar>> m_envelope;

  /// Contained surfaces of this Layer
  std::vector<std::shared_ptr<Surface>> m_surfaces = {};
};

inline const Extent& LayerBlueprint::extent() const {
  return m_extent;
}

inline const std::vector<std::shared_ptr<Surface>>& LayerBlueprint::surfaces()
    const {
  return m_surfaces;
}

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const LayerBlueprint& lbp);

}  // namespace Acts
