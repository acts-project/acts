// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"

class TGeoNode;

namespace Acts {

/// @brief ITGeoIdentierProvider
///
/// Interface class that provides an Indentifier from a TGeoNode
class ITGeoIdentifierProvider {
 public:
  /// Take a geometry context and a TGeoNode and provide an identifier
  ///
  /// @param gctx is a geometry context object
  /// @param tgnode is a TGeoNode that is translated
  virtual Identifier identify(const GeometryContext& gctx,
                              const TGeoNode& tgnode) const = 0;
};

}  // namespace Acts
