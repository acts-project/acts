// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/StaticBlueprintNode.hpp"

namespace Acts {

StaticBlueprintNode::StaticBlueprintNode(const std::string& name,
                                         std::unique_ptr<TrackingVolume> volume)
    : BlueprintNode(name), m_volume(std::move(volume)) {}

Volume& StaticBlueprintNode::build() {
  if (!m_volume) {
    throw std::runtime_error("Volume is not set");
  }

  return *m_volume;
}

void StaticBlueprintNode::connect(TrackingVolume& parent) {
  if (!m_volume) {
    throw std::runtime_error("Volume is not set");
  }

  for (auto& child : children()) {
    child->connect(*m_volume);
  }
  parent.addVolume(std::move(m_volume));
}

void StaticBlueprintNode::connect() {
  if (!m_volume) {
    throw std::runtime_error("Volume is not set");
  }

  for (auto& child : children()) {
    child->connect(*m_volume);
  }
}

}  // namespace Acts
