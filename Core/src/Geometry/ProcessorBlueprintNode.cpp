// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ProcessorBlueprintNode.hpp"

#include "Acts/Utilities/GraphViz.hpp"

#include <functional>

namespace Acts::Experimental {

namespace detail {
class ProcessorBlueprintNodeImpl {
 public:
  std::string m_name = "Unknown";

  ProcessorBlueprintNode::BuildFunction m_build;
  ProcessorBlueprintNode::ConnectFunction m_connect;
  ProcessorBlueprintNode::FinalizeFunction m_finalize;
};
}  // namespace detail

ProcessorBlueprintNode::ProcessorBlueprintNode() {
  m_impl = std::make_unique<detail::ProcessorBlueprintNodeImpl>();
}

ProcessorBlueprintNode::~ProcessorBlueprintNode() = default;

const std::string& ProcessorBlueprintNode::name() const {
  return impl().m_name;
}

ProcessorBlueprintNode& ProcessorBlueprintNode::setName(
    const std::string& name) {
  impl().m_name = name;
  return *this;
}

void ProcessorBlueprintNode::toStream(std::ostream& os) const {
  os << "ProcessorBlueprintNode(" << name() << ")";
}

Volume& ProcessorBlueprintNode::build(const BlueprintOptions& options,
                                      const GeometryContext& gctx,
                                      const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "ProcessorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "ProcessorBlueprintNode must have exactly one child");
  }

  Volume& childVolume = children().at(0).build(options, gctx, logger);
  if (impl().m_build) {
    return impl().m_build(childVolume);
  } else {
    return childVolume;
  }
}

PortalShellBase& ProcessorBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "ProcessorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "ProcessorBlueprintNode must have exactly one child");
  }

  PortalShellBase& childShell = children().at(0).connect(options, gctx, logger);

  if (impl().m_connect) {
    return impl().m_connect(options, gctx, logger, childShell);
  } else {
    return childShell;
  }
}

void ProcessorBlueprintNode::finalize(const BlueprintOptions& options,
                                      const GeometryContext& gctx,
                                      TrackingVolume& parent,
                                      const Logger& logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "ProcessorBlueprintNode must have exactly "
                           "one child, but has "
                        << children().size());
    throw std::runtime_error(
        "ProcessorBlueprintNode must have exactly one child");
  }

  children().at(0).finalize(options, gctx, parent, logger);
  if (impl().m_finalize) {
    impl().m_finalize(options, gctx, parent, logger);
  }
}

void ProcessorBlueprintNode::addToGraphviz(std::ostream& os) const {
  os << GraphViz::Node{
      .id = name(), .label = name(), .shape = GraphViz::Shape::Hexagon};
  BlueprintNode::addToGraphviz(os);
}

detail::ProcessorBlueprintNodeImpl& ProcessorBlueprintNode::impl() {
  if (!m_impl) {
    throw std::runtime_error("ProcessorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

const detail::ProcessorBlueprintNodeImpl& ProcessorBlueprintNode::impl() const {
  if (!m_impl) {
    throw std::runtime_error("ProcessorBlueprintNodeImpl is not set");
  }
  return *m_impl;
}

ProcessorBlueprintNode& ProcessorBlueprintNode::onBuild(
    ProcessorBlueprintNode::BuildFunction build) {
  impl().m_build = std::move(build);
  return *this;
}

ProcessorBlueprintNode& ProcessorBlueprintNode::onConnect(
    ProcessorBlueprintNode::ConnectFunction connect) {
  impl().m_connect = std::move(connect);
  return *this;
}

ProcessorBlueprintNode& ProcessorBlueprintNode::onFinalize(
    ProcessorBlueprintNode::FinalizeFunction finalize) {
  impl().m_finalize = std::move(finalize);
  return *this;
}

}  // namespace Acts::Experimental
