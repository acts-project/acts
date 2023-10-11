// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/ComponentBuilderProxy.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

namespace Experimental {

namespace detail {
namespace CylindricalComponentProxyHelper {

/// @brief This method helps to sort the children in
/// a container setup to according to their position/bounds
/// and the given binning
///
/// @param container [in out] the contqainer that is sorted
void sortChildren(ComponentBuilderProxy::ContainerProxy& container);

/// @brief This method is a dedicated helper method for z binning
/// @param proxy the proxy from which the z positions are extracted
/// @return a pair of positions along and opposite the z axis
std::array<Acts::Vector3, 2u> extractVolumeEndsZ(
    const ComponentBuilderProxy& proxy);

/// @brief Helper method to create a proxy for a volume in Z
///
/// @param gapName The name of the gap volume
/// @param innerR The inner radius of the gap volume
/// @param outerR The outer radious of the gap volume
/// @param transform The transform of the mother proxy
/// @param z0 the position center - half length * loca z axis
/// @param z1 the position center + half length * local z axis
/// @param logLevel the logging output level of the created builder tools
///
/// @note todo add sectoral setup
///
/// @return a newly created proxy
std::shared_ptr<ComponentBuilderProxy> createGapProxyZ(
    const std::string& gapName, ActsScalar innerR, ActsScalar outerR,
    const Transform3& transform, const Vector3& z0, const Vector3& z1,
    Logging::Level logLevel = Logging::INFO);

/// @brief Helper method to create a proxy for a volume in Z
///
/// @param gapName The name of the gap volume
/// @param innerR The inner radius of the gap volume
/// @param outerR The outer radious of the gap volume
/// @param halfZ The half length of the gap volume
/// @param transform The transform of the mother proxy
/// @param logLevel the logging output level of the created builder tools
///
/// @note todo add sectoral setup
///
/// @return a newly created proxy
std::shared_ptr<ComponentBuilderProxy> createGapProxyR(
    const std::string& gapName, ActsScalar innerR, ActsScalar outerR,
    ActsScalar halfZ, const Transform3& transform,
    Logging::Level logLevel = Logging::INFO);

/// @brief This method is a dedicated helper method to add Gap proxies
/// to a container
///
/// @param proxy the parent proxy
/// @param logLevel the logging output level of the created builder tools
///
/// fit into the bounds of the parent proxy
void addGapProxies(ComponentBuilderProxy& proxy,
                   Logging::Level logLevel = Logging::INFO);

}  // namespace CylindricalComponentProxyHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
