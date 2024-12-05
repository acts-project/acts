// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/PortalShell.hpp"

#include <variant>

namespace Acts {

class MaterialDesignatorBlueprintNode final : public BlueprintNode {
 public:
  // @TODO: This needs cuboid volume storage as well
  // @TODO: I don't love the type
  using BinningConfig = std::variant<std::vector<
      std::tuple<CylinderVolumeBounds::Face, Experimental::ProtoBinning,
                 Experimental::ProtoBinning>>>;

  MaterialDesignatorBlueprintNode(const std::string& name) : m_name(name) {}

  const std::string& name() const override;

  void toStream(std::ostream& os) const override;

  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent, const Logger& logger) override;

  void addToGraphviz(std::ostream& os) const override;

  const std::optional<BinningConfig>& binning() const;

  MaterialDesignatorBlueprintNode& setBinning(BinningConfig binning);

 private:
  std::string m_name{};

  std::optional<BinningConfig> m_binning{};
};

}  // namespace Acts
