// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"

#include <functional>

namespace Acts::Experimental {

namespace detail {
class ProcessorBlueprintNodeImpl;
}

class ProcessorBlueprintNode final : public BlueprintNode {
 public:
  using BuildFunction = std::function<Volume&(Volume&)>;

  using ConnectFunction = std::function<PortalShellBase&(
      const BlueprintOptions&, const GeometryContext&, const Logger&,
      PortalShellBase&)>;
  using FinalizeFunction =
      std::function<void(const BlueprintOptions&, const GeometryContext&,
                         TrackingVolume&, const Logger&)>;

  ProcessorBlueprintNode();
  ~ProcessorBlueprintNode() override;

  /// @copydoc BlueprintNode::name
  const std::string& name() const override;

  ProcessorBlueprintNode& setName(const std::string& name);

  /// @copydoc BlueprintNode::toStream
  void toStream(std::ostream& os) const override;

  /// @note The function is **copied** into the node! If you pass a
  ///       callable struct, make sure you don't expect side-effects on the
  ///       copied-from object.
  /// @note Calling this with an *empty* function will unset the corresponding
  ///       processing action
  ProcessorBlueprintNode& onBuild(BuildFunction build);

  /// @note The function is **copied** into the node! If you pass a
  ///       callable struct, make sure you don't expect side-effects on the
  ///       copied-from object.
  /// @note Calling this with an *empty* function will unset the corresponding
  ///       processing action
  ProcessorBlueprintNode& onConnect(ConnectFunction connect);

  /// @note The function is **copied** into the node! If you pass a
  ///       callable struct, make sure you don't expect side-effects on the
  ///       copied-from object.
  /// @note Calling this with an *empty* function will unset the corresponding
  ///       processing action
  ProcessorBlueprintNode& onFinalize(FinalizeFunction finalize);

  /// @name BlueprintNode interface
  /// @{

  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent, const Logger& logger) override;
  /// @}

 private:
  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override;

  detail::ProcessorBlueprintNodeImpl& impl();
  const detail::ProcessorBlueprintNodeImpl& impl() const;

  std::unique_ptr<detail::ProcessorBlueprintNodeImpl> m_impl;
};

}  // namespace Acts::Experimental
