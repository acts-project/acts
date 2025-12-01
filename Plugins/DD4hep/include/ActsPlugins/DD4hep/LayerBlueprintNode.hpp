// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/LayerBlueprintNode.hpp"

#include <string>

#include <DD4hep/DetElement.h>

namespace ActsPlugins::DD4hep {

class LayerBlueprintNode : public Acts::Experimental::LayerBlueprintNode {
 public:
  //   explicit LayerBlueprintNode(const dd4hep::DetElement& detElement,
  //                               std::string_view axes)
  //       : Acts::Experimental::LayerBlueprintNode(detElement.name()) {}

  using Acts::Experimental::LayerBlueprintNode::LayerBlueprintNode;
};

}  // namespace ActsPlugins::DD4hep