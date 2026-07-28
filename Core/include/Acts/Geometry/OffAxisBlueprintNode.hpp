// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/StaticBlueprintNode.hpp"

namespace Acts {

/// Blueprint node that wraps one off-axis subtree in an on-axis cylindrical
/// envelope.
///
/// The wrapped child subtree keeps its original (possibly off-axis) placement.
/// This node computes a conservative on-axis cylinder around that subtree and
/// presents this envelope to parent cylinder stacks.
class OffAxisBlueprintNode final : public StaticBlueprintNode {
 public:
  /// Construct the off-axis wrapper.
  /// @param name Name of the wrapper volume
  /// @param axisTransform Reference transform of the co-axial parent axis
  explicit OffAxisBlueprintNode(
      std::string_view name,
      Transform3 axisTransform = Transform3::Identity());

  /// @copydoc BlueprintNode::name
  const std::string& name() const override;

  /// @copydoc BlueprintNode::build
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// Set the reference axis transform used for envelope construction.
  /// @param axisTransform Transform of the reference axis frame
  /// @return Reference to this node for chaining
  OffAxisBlueprintNode& setAxisTransform(const Transform3& axisTransform);

  /// Get the reference axis transform.
  /// @return Reference to the current axis transform
  const Transform3& axisTransform() const;

 private:
  std::string m_name;
  Transform3 m_axisTransform = Transform3::Identity();
};

}  // namespace Acts
