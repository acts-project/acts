// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// A implementation of a detector element, that is constructed from a
/// JSON description of a surface. The idea behind this is that it helps
/// importing whole tracking geometries from JSON files. In some parts of
/// the codebase, the existence of a detector element associated to a surface
/// has a specific meaning (e.g., flags surfaces as sensitive).
class JsonDetectorElement : public DetectorElementBase {
 public:
  /// Constructor from JSON surface description
  /// @param jSurface JSON object describing the surface
  /// @param thickness Thickness of the detector element
  JsonDetectorElement(const nlohmann::json &jSurface, double thickness);

  /// Return mutable reference to the surface
  /// @return Mutable reference to the associated surface
  Surface &surface() override;
  /// Return const reference to the surface
  /// @return Const reference to the associated surface
  const Surface &surface() const override;

  /// Return the thickness of the detector element
  /// @return Thickness value
  double thickness() const override;

  /// Return the transform for this detector element
  /// @param gctx Geometry context (unused for this implementation)
  /// @return Transform matrix for this detector element
  const Transform3 &localToGlobalTransform(
      const GeometryContext &gctx) const override;

  bool isSensitive() const override { return true; }

 private:
  std::shared_ptr<Surface> m_surface;
  Transform3 m_transform{};
  double m_thickness{};
};

/// @}

}  // namespace Acts
