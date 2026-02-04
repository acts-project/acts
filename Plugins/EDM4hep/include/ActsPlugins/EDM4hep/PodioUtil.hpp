// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "ActsPlugins/EDM4hep/PodioDynamicColumns.hpp"

#include <limits>
#include <memory>

#include <podio/podioVersion.h>

#if podio_VERSION_MAJOR >= 1 || \
    (podio_VERSION_MAJOR == 0 && podio_VERSION_MINOR == 99)
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>
#else
#include <podio/ROOTFrameReader.h>
#include <podio/ROOTFrameWriter.h>
#endif

namespace ActsPodioEdm {
class Surface;
}

namespace podio {
class Frame;
}

namespace ActsPlugins {
/// @addtogroup edm4hep_plugin
/// @{

namespace PodioUtil {

// We want to support podio 0.16 and 1.x for now

// See https://github.com/AIDASoft/podio/pull/549
#if podio_VERSION_MAJOR >= 1 || \
    (podio_VERSION_MAJOR == 0 && podio_VERSION_MINOR == 99)
using ROOTWriter = podio::ROOTWriter;
using ROOTReader = podio::ROOTReader;
#else
using ROOTWriter = podio::ROOTFrameWriter;
using ROOTReader = podio::ROOTFrameReader;
#endif

// See https://github.com/AIDASoft/podio/pull/553
template <typename T>
decltype(auto) getDataMutable(T&& object) {
  if constexpr (podio::version::build_version.major >= 1) {
    return std::forward<T>(object).getData();
  } else {
    return std::forward<T>(object).data();
  }
}

template <typename T>
decltype(auto) getReferenceSurfaceMutable(T&& object) {
  if constexpr (podio::version::build_version.major >= 1) {
    return std::forward<T>(object).getReferenceSurface();
  } else {
    return std::forward<T>(object).referenceSurface();
  }
}

using Identifier = std::uint64_t;
constexpr Identifier kNoIdentifier = std::numeric_limits<Identifier>::max();
constexpr int kNoSurface = -1;

/// Helper interface for converting between ACTS and PODIO types
/// @ingroup edm4hep_plugin
class ConversionHelper {
 public:
  /// Convert surface to identifier
  /// @param surface The surface to convert
  /// @return Optional identifier for the surface
  virtual std::optional<Identifier> surfaceToIdentifier(
      const Acts::Surface& surface) const = 0;
  /// Convert identifier to surface
  /// @param identifier The identifier to convert
  /// @return Pointer to the surface
  virtual const Acts::Surface* identifierToSurface(
      Identifier identifier) const = 0;

  /// Convert source link to identifier
  /// @param sl The source link to convert
  /// @return Identifier for the source link
  virtual Identifier sourceLinkToIdentifier(const Acts::SourceLink& sl) = 0;
  /// Convert identifier to source link
  /// @param identifier The identifier to convert
  /// @return Source link for the identifier
  virtual Acts::SourceLink identifierToSourceLink(
      Identifier identifier) const = 0;
};

std::shared_ptr<const Acts::Surface> convertSurfaceFromPodio(
    const ConversionHelper& helper, const ActsPodioEdm::Surface& surface);

ActsPodioEdm::Surface convertSurfaceToPodio(const ConversionHelper& helper,
                                            const Acts::Surface& surface);
}  // namespace PodioUtil

/// @}

namespace podio_detail {
/// This is used by both the track and track state container, so the
/// implementation is shared here
void recoverDynamicColumns(
    const podio::Frame& frame, const std::string& stem,
    std::unordered_map<Acts::HashedString,
                       std::unique_ptr<podio_detail::ConstDynamicColumnBase>>&
        dynamic);
}  // namespace podio_detail

}  // namespace ActsPlugins
