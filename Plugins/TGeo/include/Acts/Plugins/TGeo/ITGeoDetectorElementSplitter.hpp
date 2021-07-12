// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"

class TGeoNode;

namespace Acts {

class TGeoDetectorElement;

/// @brief ITGeoElementSplitter
///
/// Interface class that allows to define splitting of TGeoElements into
/// sub-elements
class ITGeoDetectorElementSplitter {
 public:
  /// Take a geometry context and TGeoElement and split it into sub elements
  ///
  /// @param gctx is a geometry context object
  /// @param tgnode is a TGeoNode that is translated
  ///
  /// @note If no split is performed the unsplit detector element is returned 
  virtual std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> split(
      const GeometryContext& gctx, std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const = 0;
};

}  // namespace Acts