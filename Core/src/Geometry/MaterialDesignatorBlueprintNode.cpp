// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"

#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

const std::string& MaterialDesignatorBlueprintNode::name() const {
  return m_name;
}

void MaterialDesignatorBlueprintNode::toStream(std::ostream& os) const {
  os << "MaterialDesignatorBlueprintNode(" << name() << ")";
}

Volume& MaterialDesignatorBlueprintNode::build(const Options& options,
                                               const GeometryContext& gctx,
                                               const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "MaterialDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }
  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& MaterialDesignatorBlueprintNode::connect(
    const Options& options, const GeometryContext& gctx, const Logger& logger) {
  ACTS_DEBUG(prefix() << "MaterialDesignatorBlueprintNode::connect");
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "MaterialDesignatorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }
  if (!m_binning) {
    ACTS_ERROR(prefix() << "Binning is not set");
    throw std::runtime_error("Binning is not set");
  }

  auto& shell = children().at(0).connect(options, gctx, logger);

  ACTS_DEBUG(prefix() << "Received shell from child "
                      << children().at(0).name());

  if (auto* cylShell = dynamic_cast<CylinderPortalShell*>(&shell)) {
    ACTS_DEBUG(prefix() << "Connecting cylinder shell");

    std::visit(
        overloaded{
            [&](const std::vector<std::tuple<
                    CylinderPortalShell::Face, Experimental::ProtoBinning,
                    Experimental::ProtoBinning>>& binning) {
              ACTS_DEBUG(prefix() << "Binning is set to compatible type");
              using enum CylinderVolumeBounds::Face;

              for (auto& [face, loc0, loc1] : binning) {
                if (face == OuterCylinder || face == InnerCylinder) {
                  if (loc0.binValue != BinningValue::binRPhi) {
                    ACTS_ERROR(prefix() << "Binning is not in RPhi");
                    throw std::runtime_error("Binning is not in RPhi");
                  }

                  if (loc1.binValue != BinningValue::binZ) {
                    ACTS_ERROR(prefix() << "Binning is not in Z");
                    throw std::runtime_error("Binning is not in Z");
                  }
                }

                if (face == PositiveDisc || face == NegativeDisc) {
                  if (loc0.binValue != BinningValue::binPhi) {
                    ACTS_ERROR(prefix() << "Binning is not in Phi");
                    throw std::runtime_error("Binning is not in Phi");
                  }
                  if (loc1.binValue != BinningValue::binR) {
                    ACTS_ERROR(prefix() << "Binning is not in R");
                    throw std::runtime_error("Binning is not in R");
                  }
                }

                Experimental::BinningDescription desc{.binning = {loc0, loc1}};
                ACTS_DEBUG(prefix() << "~> Assigning proto binning "
                                    << desc.toString() << " to face " << face);

                auto material =
                    std::make_shared<ProtoGridSurfaceMaterial>(std::move(desc));

                auto portal = cylShell->portal(face);
                if (portal == nullptr) {
                  ACTS_ERROR(prefix() << "Portal is nullptr");
                  throw std::runtime_error("Portal is nullptr");
                }
                portal->surface().assignSurfaceMaterial(std::move(material));
              }
            },
            [&](const auto& /*binning*/) {
              ACTS_ERROR(prefix() << "Binning is set to unknown type");
              throw std::runtime_error("Unknown binning type");
            },
        },
        m_binning.value());
  }
  // @TODO: Handle cuboid volume shell
  else {
    ACTS_ERROR(prefix() << "Shell is not supported");
    throw std::runtime_error("Shell is not supported");
  }

  return shell;
}

void MaterialDesignatorBlueprintNode::finalize(const Options& options,
                                               TrackingVolume& parent,
                                               const Logger& logger) {
  if (children().size() != 1) {
    throw std::runtime_error(
        "MaterialDesignatorBlueprintNode must have exactly one child");
  }
  return children().at(0).finalize(options, parent, logger);
}

void MaterialDesignatorBlueprintNode::addToGraphviz(std::ostream& os) const {
  os << GraphViz::Node{.id = name(), .shape = GraphViz::Shape::InvTrapezium};
  BlueprintNode::addToGraphviz(os);
}

}  // namespace Acts
