// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoBoxConverter.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoTrdConverter.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoTubeConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>
#include <Acts/Surfaces/LineBounds.hpp>

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeShift.h>

namespace Acts {

namespace detail {

struct GeoShiftConverter {
  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoBox The GeoBox to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  template <typename ContainedShape, typename Converter, typename Surface, typename Bounds>
  std::tuple<std::shared_ptr<GeoModelDetectorElement>, std::shared_ptr<Acts::Surface>>
  impl(const GeoFullPhysVol& geoFPV, const GeoShapeShift& geoShift,
       const Transform3& absTransform, bool sensitive) const {
    auto trd = dynamic_cast<const ContainedShape*>(geoShift.getOp());

    if (trd == nullptr) {
      return {nullptr, nullptr};
    }

    const Transform3& shift = geoShift.getX();
    auto [el, surface] = Converter{}(geoFPV, *trd, absTransform * shift, sensitive);

    // Use knowledge from GeoTrdConverter to make shared bounds object
    const auto& bounds = static_cast<const Bounds&>(surface->bounds());
    auto sharedBounds = std::make_shared<const Bounds>(bounds);
    // std::cout << "     Extracted bounds params: " << sharedBounds->values()[0] << ", " << sharedBounds->values()[1] << ", " << sharedBounds->values()[2] << std::endl;

    // TODO this procedure could be stripped from all converters because it is
    // pretty generic
    if (!sensitive) {
      auto newSurface = Surface::template makeShared<Surface>(
          surface->transform({}), sharedBounds);
      return std::make_tuple(nullptr, newSurface);
    }

    auto newEl = GeoModelDetectorElement::createDetectorElement<Surface>(
        el->physicalVolume(), sharedBounds, el->transform({}),
        el->thickness());
    auto newSurface = newEl->surface().getSharedPtr();
    return std::make_tuple(newEl, newSurface);
  }

  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoBox The GeoBox to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  std::tuple<std::shared_ptr<GeoModelDetectorElement>, std::shared_ptr<Surface>>
  operator()(const GeoFullPhysVol& geoFPV, const GeoShapeShift& geoShift,
             const Transform3& absTransform, bool sensitive) const {
    std::tuple<std::shared_ptr<GeoModelDetectorElement>,
               std::shared_ptr<Surface>>
        r;
    auto& [el, surface] = r;

    r = impl<GeoTrd, detail::GeoTrdConverter, PlaneSurface, TrapezoidBounds>(
        geoFPV, geoShift, absTransform, sensitive);

    if (surface) {
      // Acts::GeometryContext gctx{};
      // std::cout << "Final surface: " << std::tie(*surface, gctx) << std::endl;
      return r;
    }

    r = impl<GeoBox, detail::GeoBoxConverter, PlaneSurface, RectangleBounds>(
        geoFPV, geoShift, absTransform, sensitive);

    if (surface) {
      return r;
    }

    // For now this does straw by default
    r = impl<GeoTube, detail::GeoTubeConverter, StrawSurface, LineBounds>(
        geoFPV, geoShift, absTransform, sensitive);

    if (surface) {
      return r;
    }

    return {nullptr, nullptr};
  }
};
}  // namespace detail

/// @brief The GeoShift + Trd converter
///
/// This is a dedicated converter for GeoBox shapes
using GeoShiftConverter =
    detail::GenericGeoShapeConverter<GeoShapeShift, detail::GeoShiftConverter>;

}  // namespace Acts
