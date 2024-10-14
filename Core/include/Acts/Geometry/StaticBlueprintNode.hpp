// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

class StaticBlueprintNode : public BlueprintNode {
 public:
  StaticBlueprintNode(std::unique_ptr<TrackingVolume> volume);

  const std::string& name() const override;

  Volume& build(const Options& options, const GeometryContext& gctx,

                const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const Options& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const Options& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  virtual StaticBlueprintNode& setNavigationPolicyFactory(
      std::shared_ptr<NavigationPolicyFactory> navigationPolicyFactory);

  const NavigationPolicyFactory* navigationPolicyFactory() const;

 protected:
  void addToGraphviz(std::ostream& os) const override;

  std::unique_ptr<TrackingVolume> m_volume;

  std::unique_ptr<PortalShellBase> m_shell;

  std::shared_ptr<NavigationPolicyFactory> m_navigationPolicyFactory = nullptr;
};

}  // namespace Acts
