// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Experimental/SurfaceLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

class Surface;
class VolumeBounds;

using InternalSurfaceLinks =
    std::tuple<SurfaceLinks, std::vector<SurfaceLinks>>;

using SurfaceLinksGenerator =
    std::function <
    InternalSurfaceLinks(const GeometryContext&,
                         const std::vector<std::shared_ptr<Surface>>&,
                         const VolumeBounds*)>;

/// A non-specified surface link generator that creates all dummy surface links
/// to all internal surfaces
struct AllInternalSurfaces {
  /// The number of boundaries to be filled
  size_t boundaries;

  /// This is a standard surface link creator that provides all surfaces
  /// without local search
  ///
  InternalSurfaceLinks operator()(
      const GeometryContext& /*ignored*/,
      const std::vector<std::shared_ptr<Surface>>& /*ignored*/,
      const VolumeBounds* /*ignored*/) const {
    return {AllSurfaces{},
            std::vector<SurfaceLinks>(boundaries, AllSurfaces{})};
  }
};

/// The blueprint for internal structure
///
/// This class is designed to describe how layers are built
class InternalBlueprint {
 public:
  /// (Slow) Constructor with all/part of the surfaces
  ///
  /// @param gctx The geometry context
  /// @param surfaces The surfaces in to be contained by this layer volume
  /// @param surfaceLinksGenerator The surface links gnerator
  /// @param binValueEnvelopes The list of extent binning values, empty means all
  /// @param name The name of the blueprint
  InternalBlueprint(const GeometryContext& gctx,
                    const std::vector<std::shared_ptr<Surface>>& surfaces = {},
                    const SurfaceLinksGenerator& surfaceLinksGenerator =
                        AllInternalSurfaces{},
                    const std::vector<std::pair<BinningValue, Envelope>>&
                        binValueEnvelopes = {},
                    const std::string& name = "Unnamed");

  /// (Fast) Constructor with all/part of the surfaces
  ///
  /// @note this is meant for cases where the extents have been precomputed
  ///
  /// @param surfacesExtents The surfaces in to be contained by this layer volume
  /// @param surfaceLinksGenerator The surface links gnerator
  /// @param name The name of the blueprint
  InternalBlueprint(
      const std::vector<std::pair<std::shared_ptr<Surface>, GeometricExtent>>&
          surfacesExtent = {},
      const SurfaceLinksGenerator& surfaceLinksGenerator =
          AllInternalSurfaces{},
      const std::string& name = "Unnamed");

  InternalBlueprint() = default;

  /// @return the current extent
  GeometricExtent& extent();

  /// @return the current extent - non-const access
  const GeometricExtent& extent() const;

  /// @return the vector of surfaces
  const std::vector<std::shared_ptr<Surface>>& surfaces() const;

  /// @return the surface links object
  const SurfaceLinksGenerator& surfaceLinksGenerator() const;

  /// @return the name of the layer
  const std::string& name() const;

  /// Output to ostream
  /// @param sl the input ostream
  std::ostream& toStream(std::ostream& sl) const;

 private:
  /// Contained surfaces of this Layer
  std::vector<std::shared_ptr<Surface>> m_surfaces = {};

  /// The surface links object
  SurfaceLinksGenerator m_surfaceLinksGenerator = AllInternalSurfaces{};

  /// The actual extent (measured from surfaces)
  GeometricExtent m_extent = GeometricExtent{};

  /// The name of the blue print
  std::string m_name;
};

inline GeometricExtent& InternalBlueprint::extent() {
  return m_extent;
}

inline const GeometricExtent& InternalBlueprint::extent() const {
  return m_extent;
}

inline const std::vector<std::shared_ptr<Surface>>&
InternalBlueprint::surfaces() const {
  return m_surfaces;
}

inline const SurfaceLinksGenerator& InternalBlueprint::surfaceLinksGenerator()
    const {
  return m_surfaceLinksGenerator;
}

inline const std::string& InternalBlueprint::name() const {
  return m_name;
}

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const InternalBlueprint& lbp);

}  // namespace Acts
