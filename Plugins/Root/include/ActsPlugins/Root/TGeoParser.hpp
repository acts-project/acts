// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TGeoMatrix.h"

class TGeoNode;
class TGeoVolume;

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// @brief TGeoParser is a helper struct that
/// walks recursively through a TGeometry and selects by
/// string comparison the TGeoNodes that match the criteria
///
/// It also buils up the global transform for the conversion
/// into an ACTS Surface
struct TGeoParser {
  /// Type alias for parsing range as min/max bounds pair
  using ParseRange = std::pair<double, double>;

  /// Structure holding a selected TGeo node with its global transform
  struct SelectedNode {
    /// The selected geo node
    const TGeoNode* node = nullptr;
    /// The transform to global
    std::unique_ptr<TGeoMatrix> transform = nullptr;
  };

  /// @brief Nested state struct
  ///
  /// This is needed for the recursive parsing of the
  /// geometry, it collects the information during the parsing process
  /// and keeps track of the built up transform
  struct State {
    /// The current ROOT geometry volume being parsed
    TGeoVolume* volume = nullptr;
    /// The current ROOT geometry node being parsed
    TGeoNode* node = nullptr;
    /// Flag indicating if parsing is on the selected geometry branch
    bool onBranch = false;
    /// Collection of nodes selected during parsing
    std::vector<SelectedNode> selectedNodes = {};
  };

  /// @brief Nested configuration struct
  ///
  /// This contains the parsing configuration
  struct Options {
    /// Identify the volume by name
    std::vector<std::string> volumeNames = {};
    /// Identify the sensor(s)/target(s) by name
    std::vector<std::string> targetNames = {};
    /// The local axis definition of TGeo object wrt Acts::Surface
    std::string localAxes = "XYZ";
    /// Scaling from TGeo to ROOT
    double unit = 1 * Acts::UnitConstants::cm;
    /// Parse restrictions, several can apply
    std::vector<std::pair<Acts::AxisDirection, ParseRange> > parseRanges = {};
  };

  /// The parsing module, it takes the top Volume and recursively steps down
  /// @param state [out] The parsing state configuration, passed through
  /// @param options [in] The parsing options as required
  /// @param gmatrix The current built-up transform to global at this depth
  static void select(State& state, const Options& options,
                     const TGeoMatrix& gmatrix = TGeoIdentity("ID"));

  /// Simple utility function that recursively finds the node by the volume name
  /// in the tgeo branch.
  /// @param currentNode [in] the pointer to the current node in the branch
  /// @param volumeName  [in] the name of the volume to be searched for
  /// @return the pointer to the node corresponding to the volume
  static TGeoNode* findNodeRecursive(TGeoNode* currentNode,
                                     const char* volumeName);
};
/// @}
}  // namespace ActsPlugins
