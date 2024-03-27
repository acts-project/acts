// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/BlueprintNode.hpp"

#include <ostream>

namespace Acts {

void BlueprintNode::toStream(std::ostream& os) const {
  os << "BlueprintNode( " << name() << ")";
}

void BlueprintNode::addChild(std::unique_ptr<BlueprintNode> child) {
  child->m_depth = m_depth + 1;
  m_children.push_back(std::move(child));
}

BlueprintNode::MutableChildRange BlueprintNode::children() {
  return MutableChildRange{m_children};
}

BlueprintNode::ChildRange BlueprintNode::children() const {
  return ChildRange{m_children};
}

std::size_t BlueprintNode::depth() const {
  return m_depth;
}

std::string BlueprintNode::prefix() const {
  return std::string(m_depth * 2, ' ') + "[" + name() + "]: ";
}

void BlueprintNode::visualize(IVisualization3D& vis,
                              const GeometryContext& gctx) const {
  for (const auto& child : children()) {
    child.visualize(vis, gctx);
  }
}

}  // namespace Acts
