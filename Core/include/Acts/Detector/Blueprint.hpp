// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {
namespace Experimental {

class IGeometryIdGenerator;
class IInternalStructureBuilder;
class IRootVolumeFinderBuilder;

/// A Blueprint is an instruction tree that allows you to defina a tree sequence
/// of volume building using the provided tools.
///
/// It follows tree nomenclature and can define:
///
/// - a root node (the top of the tree)
/// - a branch node (also called inner node)
/// - leaf node (also called terminal node)
///
/// Leaf nodes can have internal builders attached, while all the external
/// builders will be created when the blueprint is interpreted.
namespace Blueprint {

struct Node final {
  /// Branch constructor
  ///
  /// @param n name of the node
  /// @param t the transform
  /// @param bt the boundary type
  /// @param bv the boundary values
  /// @param bss the binning values
  /// @param cs the children of the node
  Node(const std::string& n, const Transform3& t, VolumeBounds::BoundsType bt,
       const std::vector<ActsScalar>& bv, const std::vector<BinningValue>& bss,
       std::vector<std::unique_ptr<Node>> cs = {})
      : name(n),
        transform(t),
        boundsType(bt),
        boundaryValues(bv),
        children(std::move(cs)),
        binning(bss) {
    for_each(children.begin(), children.end(),
             [this](std::unique_ptr<Node>& c) { c->parent = this; });
  }

  /// Leaf constructor
  ///
  /// @param n name of the node
  /// @param t the transform
  /// @param bt the boundary type
  /// @param bv the boundary values
  /// @param isb the internal structure builder (optional)
  Node(const std::string& n, const Transform3& t, VolumeBounds::BoundsType bt,
       const std::vector<ActsScalar>& bv,
       std::shared_ptr<const IInternalStructureBuilder> isb = nullptr)
      : name(n),
        transform(t),
        boundsType(bt),
        boundaryValues(bv),
        internalsBuilder(std::move(isb)) {}

  /// Name identification of this node
  std::string name = "";
  /// Transform definition of this node
  Transform3 transform = Transform3::Identity();
  /// Boundary definition of this node
  VolumeBounds::BoundsType boundsType = VolumeBounds::eOther;
  /// The boundary type
  std::vector<ActsScalar> boundaryValues = {};
  /// Parent node - nullptr for root only
  const Node* parent = nullptr;
  /// Branch definitions: children
  std::vector<std::unique_ptr<Node>> children = {};
  /// Branch definition binning
  std::vector<BinningValue> binning = {};

  /// Builders and helper tools that can be attached
  std::shared_ptr<const IRootVolumeFinderBuilder> rootVolumeFinderBuilder =
      nullptr;
  /// Geometry id generator
  std::shared_ptr<const IGeometryIdGenerator> geoIdGenerator = nullptr;
  /// Internal structure builder - for leaf nodes
  std::shared_ptr<const IInternalStructureBuilder> internalsBuilder = nullptr;

  /// @brief Check if it is a leaf node
  bool isLeaf() const { return children.empty(); }

  /// @brief Check is it is a root
  bool isRoot() const { return parent == nullptr; }

  /// @brief Method to add a child to this branch
  /// @param c the child to be added
  void add(std::unique_ptr<Node> c) {
    c->parent = this;
    children.push_back(std::move(c));
  }

  /// @brief  Turn into a dot output
  template <typename stream_type>
  void dotStream(stream_type& ss,
                 const std::string& graphName = "blueprint") const {
    if (isRoot()) {
      ss << "digraph " << graphName << " {" << '\n';
      ss << name
         << " [shape=\"circle\";style=\"filled\";fillcolor=\"darkorange\"];"
         << '\n';

    } else if (isLeaf()) {
      std::string color =
          (internalsBuilder != nullptr) ? "darkolivegreen1" : "darkolivegreen3";

      ss << name << " [shape=\"box\";style=\"filled\";fillcolor=\"";
      ss << color << "\"];" << '\n';
    } else {
      ss << name << " [shape=\"diamond\"];" << '\n';
    }

    ss << name << " [label=\"" << name << "\"];" << '\n';
    for (const auto& c : children) {
      ss << name << " -> " << c->name << ";" << '\n';
      c->dotStream(ss);
    }
    if (children.empty()) {
      ss << name + "_shape"
         << " [shape=\"cylinder\";style=\"filled\";fillcolor=\"lightgrey\"];"
         << '\n';
      ss << name << " -> " << name + "_shape"
         << ";" << '\n';
    }

    if (internalsBuilder != nullptr) {
      ss << name + "_int"
         << " [shape=\"doubleoctagon\";style=\"filled\";fillcolor="
            "\"cadetblue1\"];"
         << '\n';
      ss << name << " -> " << name + "_int"
         << ";" << '\n';
    }

    if (geoIdGenerator != nullptr) {
      ss << name + "_geoID"
         << " [shape=\"note\";style=\"filled\";fillcolor=\"azure\"];" << '\n';
      ss << name << " -> " << name + "_geoID"
         << ";" << '\n';
    }

    if (rootVolumeFinderBuilder != nullptr) {
      ss << name + "_roots"
         << " [shape=\"Msquare\";style=\"filled\";fillcolor=\"darkkhaki\"];"
         << '\n';
      ss << name << " -> " << name + "_roots"
         << ";" << '\n';
    }

    if (isRoot()) {
      ss << "}" << '\n';
    }
  }
};

}  // namespace Blueprint
}  // namespace Experimental
}  // namespace Acts
