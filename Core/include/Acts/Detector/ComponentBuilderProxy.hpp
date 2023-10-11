// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <memory>
#include <string>

namespace Acts {
namespace Experimental {

class IDetectorComponentBuilder;

/// @brief A container proxy struct that can be used to build up
/// a tree graph of detector components
class ComponentBuilderProxy
    : public std::enable_shared_from_this<ComponentBuilderProxy> {
 protected:
  /// Createing a ComponentBuilderProxy is not allowed
  /// @note use : createRootProxy(), createChildVolumeProxy() or
  /// createChildContainerProxy()
  ComponentBuilderProxy() = default;

 public:
  /// @brief Private constructor to be used for the proxy creation
  /// @param n The name of the proxy
  ComponentBuilderProxy(const std::string& n) : name(n){};

  /// @brief  Definition of a volume
  struct VolumeProxy {
    /// If fully defined, the builder is set
    std::shared_ptr<const IDetectorComponentBuilder> builder = nullptr;
  };

  /// @brief  Definition of a container
  struct ContainerProxy {
    /// The binning instruction
    std::vector<BinningValue> binning = {};
    /// The contained component builders, in case of a container builder
    std::vector<std::shared_ptr<ComponentBuilderProxy>> children = {};
  };

  /// The proxy name - this needs to be uniuque
  std::string name = "";

  /// The Parent proxy, optionally set
  std::shared_ptr<ComponentBuilderProxy> parent = nullptr;

  /// The volume proxy, optionally set
  std::variant<VolumeProxy, ContainerProxy> holder = VolumeProxy{};

  /// The boundary type
  VolumeBounds::BoundsType boundsType = VolumeBounds::eOther;

  /// The boundary values - to be interpreted together with boundary type
  std::vector<ActsScalar> boundaryValues = {};

  /// The transform
  Transform3 transform = Transform3::Identity();

  /// Request the builder of this object
  std::shared_ptr<const IDetectorComponentBuilder> builder() const {
    if (std::holds_alternative<VolumeProxy>(this->holder)) {
      return std::get<VolumeProxy>(this->holder).builder;
    } else {
      std::runtime_error(
          "ComponentBuilderProxy::builder() - "
          "Cannot return a builder for a container proxy yet, "
          "the detecor setup is not complete.");
    }
    return nullptr;
  }

  /// @brief Create the Root volume proxy
  ///
  /// @param n The name of the Root proxy
  /// @param t The transform of the Root proxy
  /// @param bt The boundary type of the Root proxy
  /// @param bv The boundary values of the Root proxy
  /// @param b The builder of the Root proxy - in this case
  /// it is a single volume builder case
  ///
  /// @return The created Root proxy as a shared_ptr
  static std::shared_ptr<ComponentBuilderProxy> createRootProxy(
      std::string n, Transform3 t = Transform3::Identity(),
      VolumeBounds::BoundsType bt = VolumeBounds::eOther,
      const std::vector<ActsScalar>& bv = {},
      std::vector<BinningValue> binning = {},
      std::shared_ptr<const IDetectorComponentBuilder> b = nullptr) {
    auto proxy = std::make_shared<ComponentBuilderProxy>(n);
    proxy->transform = t;
    proxy->boundsType = bt;
    proxy->boundaryValues = bv;
    if (b != nullptr) {
      proxy->holder = VolumeProxy{b};
    } else {
      proxy->holder = ContainerProxy{binning};
    }
    return proxy;
  }

  /// @brief Create a child volume proxy
  ///
  /// @param n The name of the child volume proxy
  /// @param t The transform of the child volume proxy
  /// @param bt The boundary type of the child volume proxy
  /// @param bv The boundary values of the child volume proxy
  /// @param b The builder of the child volume proxy
  std::shared_ptr<ComponentBuilderProxy> addChildVolumeProxy(
      std::string n, Transform3 t = Transform3::Identity(),
      VolumeBounds::BoundsType bt = VolumeBounds::eOther,
      const std::vector<ActsScalar>& bv = {},
      std::shared_ptr<const IDetectorComponentBuilder> b = nullptr) {
    auto proxy = std::make_shared<ComponentBuilderProxy>(n);
    proxy->transform = t;
    proxy->boundsType = bt;
    proxy->boundaryValues = bv;
    proxy->holder = VolumeProxy{b};
    if (std::holds_alternative<ContainerProxy>(this->holder)) {
      std::get<ContainerProxy>(this->holder).children.push_back(proxy);
    } else {
      throw std::runtime_error(
          "ComponentBuilderProxy::addChildVolumeProxy() - "
          "Cannot add a child volume proxy to a non-container proxy");
    }
    proxy->parent = this->shared_from_this();
    return proxy;
  }

  /// @brief Create a child container proxy
  ///
  /// @param n The name of the child container proxy
  /// @param t The transform of the child container proxy
  /// @param bt The boundary type of the child container proxy
  /// @param bv The boundary values of the child container proxy
  std::shared_ptr<ComponentBuilderProxy> addChildContainerProxy(
      std::string n, Transform3 t = Transform3::Identity(),
      VolumeBounds::BoundsType bt = VolumeBounds::eOther,
      const std::vector<ActsScalar>& bv = {},
      const std::vector<BinningValue>& binning = {}) {
    auto proxy = std::make_shared<ComponentBuilderProxy>(n);
    proxy->transform = t;
    proxy->boundsType = bt;
    proxy->boundaryValues = bv;
    proxy->holder = ContainerProxy{binning};
    if (std::holds_alternative<ContainerProxy>(this->holder)) {
      std::get<ContainerProxy>(this->holder).children.push_back(proxy);
    } else {
      throw std::runtime_error(
          "ComponentBuilderProxy::addChildVolumeProxy() - "
          "Cannot add a child container proxy to a non-container proxy");
    }
    proxy->parent = this->shared_from_this();
    return proxy;
  }
};

}  // namespace Experimental
}  // namespace Acts
