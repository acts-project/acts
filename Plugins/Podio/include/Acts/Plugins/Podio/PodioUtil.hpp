// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Podio/PodioDynamicColumns.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <limits>
#include <memory>

#include <podio/podioVersion.h>
#if PODIO_VERSION_MAJOR >= 1
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

namespace Acts {
namespace PodioUtil {

// We want to support podio 0.16 and 1.x for now

// See https://github.com/AIDASoft/podio/pull/549
#if PODIO_VERSION_MAJOR >= 1
using ROOTWriter = podio::ROOTWriter;
using ROOTReader = podio::ROOTReader;

template <typename T>
auto getData(T& object) {
  return object.getData();
}

#else
using ROOTWriter = podio::ROOTWriter;
using ROOTReader = podio::ROOTReader;

#endif

// See https://github.com/AIDASoft/podio/pull/553
template <typename T>
decltype(auto) getData(const T& object) {
  return object.getData();
}

template <typename T>
decltype(auto) getDataMutable(T& object) {
#if PODIO_VERSION_MAJOR >= 1
  return object.getData();
#else
  return object.data();
#endif
}

using Identifier = std::uint64_t;
constexpr Identifier kNoIdentifier = std::numeric_limits<Identifier>::max();
constexpr int kNoSurface = -1;

// @TODO: We might want to consider making this a type erased type that's not an interface
class ConversionHelper {
 public:
  virtual std::optional<Identifier> surfaceToIdentifier(
      const Surface& surface) const = 0;
  virtual const Surface* identifierToSurface(Identifier identifier) const = 0;

  virtual Identifier sourceLinkToIdentifier(const SourceLink& sl) = 0;
  virtual SourceLink identifierToSourceLink(Identifier identifier) const = 0;
};

std::shared_ptr<const Surface> convertSurfaceFromPodio(
    const ConversionHelper& helper, const ActsPodioEdm::Surface& surface);

ActsPodioEdm::Surface convertSurfaceToPodio(const ConversionHelper& helper,
                                            const Acts::Surface& surface);
}  // namespace PodioUtil

namespace podio_detail {
/// This is used by both the track and track state container, so the
/// implementation is shared here
void recoverDynamicColumns(
    const podio::Frame& frame, const std::string& stem,
    std::unordered_map<HashedString,
                       std::unique_ptr<podio_detail::ConstDynamicColumnBase>>&
        dynamic);
}  // namespace podio_detail
}  // namespace Acts
