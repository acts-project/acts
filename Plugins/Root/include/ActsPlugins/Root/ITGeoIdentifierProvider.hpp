// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsPlugins/Root/TGeoDetectorElement.hpp"

class TGeoNode;

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// @brief ITGeoIdentierProvider
///
/// Interface class that provides an Identifier from a TGeoNode
class ITGeoIdentifierProvider {
 public:
  /// Virtual destructor
  virtual ~ITGeoIdentifierProvider() = default;

  /// Take a geometry context and a TGeoNode and provide an identifier
  ///
  /// @param gctx is a geometry context object
  /// @param tgnode is a TGeoNode that is translated
  /// @return The detector element identifier for the given TGeo node
  virtual TGeoDetectorElement::Identifier identify(
      const Acts::GeometryContext& gctx, const TGeoNode& tgnode) const = 0;
};

/// @}
}  // namespace ActsPlugins
