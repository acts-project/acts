// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts::Experimental {

/// Wraps a single child blueprint node and pads it into a larger volume whose
/// dimensions are evaluated at construction time using the child's extent plus
/// a configurable envelope.
/// @note This node can only have a single child. This is not an error during
///       tree building, but during geometry construction.
/// It defers most of the functionality to @ref Acts::Experimental::StaticBlueprintNode,
/// and only implements the build phase to perform the padding.
class PadBlueprintNode final : public StaticBlueprintNode {
 public:
  /// Main constructor for the padding node.
  /// @param name The name of the padded volume.
  /// @param envelope The envelope to apply to the child node's extent to create the padded volume.
  explicit PadBlueprintNode(const std::string& name,
                            const ExtentEnvelope& envelope);

  ~PadBlueprintNode() override = default;

  /// @copydoc BlueprintNode::build
  /// Build-phase of the blueprint construction. Returns the padded volume.
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// Get the tracking volume for this node.
  /// @return Pointer to the tracking volume (may be nullptr)
  TrackingVolume* trackingVolume() const;

 private:
  ExtentEnvelope m_envelope;
  std::string m_name;
};

}  // namespace Acts::Experimental
