// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"

namespace Acts {

class MaterialDesignatorBlueprintNode final : public BlueprintNode {
 public:
  const std::string& name() const override {
    const static std::string s_name = "Material";
    return s_name;
  }

  void toStream(std::ostream& os) const override {
    os << "MaterialDesignatorBlueprintNode(" << name() << ")";
  }

  Volume& build(const Options& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override {
    if (children().size() != 1) {
      throw std::runtime_error(
          "MaterialDesignatorBlueprintNode must have exactly one child");
    }
    return children().at(0).build(options, gctx, logger);
  }

  PortalShellBase& connect(
      const Options& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override {
    if (children().size() != 1) {
      throw std::runtime_error(
          "MaterialDesignatorBlueprintNode must have exactly one child");
    }
    return children().at(0).connect(options, gctx, logger);
  }

  void finalize(const Options& options, TrackingVolume& parent,
                const Logger& logger) override {
    if (children().size() != 1) {
      throw std::runtime_error(
          "MaterialDesignatorBlueprintNode must have exactly one child");
    }
    return children().at(0).finalize(options, parent, logger);
  }

  void addToGraphviz(std::ostream& os) const override {
    os << GraphViz::Node{.id = name(), .shape = GraphViz::Shape::InvTrapezium};
    BlueprintNode::addToGraphviz(os);
  }

 private:
  std::string m_name;
  std::string m_material;
};

}  // namespace Acts
