// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPlugins/Detray/DetrayMetadata.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include <detray/core/detector.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsPlugins {

/// @ingroup detray_plugin
/// @brief Converter from an ACTS TrackingGeometry to a detray detector
///
/// The geometry conversion is a two-step process: first the ACTS geometry is
/// converted into detray payloads by a @c DetrayPayloadConverter, then those
/// payloads are used to build the actual detray detector. This class drives
/// the second step and lets the call site fully customize the first step by
/// supplying its own configured @c DetrayPayloadConverter instance through the
/// @c Config.
///
/// @c convert is a member template over the detray metadata type. It is
/// explicitly instantiated (and declared `extern template`) for the closed set
/// of supported metadata types listed in @ref DetrayMetadata.hpp, so that the
/// heavy detray detector-building code is compiled only once per metadata in a
/// dedicated translation unit. Experiment code may still convert to a custom
/// metadata by including @c DetrayGeometryConverter.ipp and instantiating
/// @c convert with the desired metadata type.
class DetrayGeometryConverter {
 public:
  /// @brief Configuration for the geometry converter
  struct Config {
    /// The payload converter used to turn the ACTS geometry into detray
    /// payloads. Supplying it here allows the call site to fully customize the
    /// payload conversion (e.g. the beampipe volume, sensitive surface
    /// strategy or the navigation/material dispatchers).
    std::shared_ptr<const DetrayPayloadConverter> payloadConverter;

    /// Whether to convert material information from ACTS to detray
    bool convertMaterial = true;

    /// Whether to convert surface grid information from ACTS to detray
    bool convertSurfaceGrids = true;
  };

  /// @brief Combined result of a geometry conversion
  ///
  /// Bundles the built detray detector together with the detray volume/surface
  /// name map.
  /// @tparam metadata_t the detector metadata type
  template <typename metadata_t>
  struct DetrayGeometry {
    /// The built detray detector
    std::shared_ptr<detray::detector<metadata_t>> detector;

    /// The detray volume and surface name map
    detray::name_map names;
  };

  /// Constructor
  /// @param config Configuration object
  /// @param logger Logger instance
  explicit DetrayGeometryConverter(
      Config config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DetrayGeometryConverter", Acts::Logging::INFO))
      : m_cfg(std::move(config)), m_logger(std::move(logger)) {
    if (m_cfg.payloadConverter == nullptr) {
      throw std::invalid_argument(
          "DetrayGeometryConverter: payloadConverter must be set");
    }
  }

  /// @brief Convert an ACTS TrackingGeometry into a detray detector
  /// @tparam metadata_t the detector metadata type to build
  ///
  /// @param mr the memory resource to use for the detray detector construction
  /// @param gctx the geometry context
  /// @param trackingGeometry the ACTS tracking geometry to convert
  /// @param detectorName the name to set for the detray detector (optional, if
  ///     not set, it will be taken from the payloads or defaulted to empty)
  ///
  /// This method performs the following steps:
  /// 1. It converts the ACTS tracking geometry into detray payloads using the
  ///    configured DetrayPayloadConverter.
  /// 2. It builds a detray detector from the converted payloads using the
  ///    detray::detector_builder.
  ///
  /// @note The definition lives in @c DetrayGeometryConverter.ipp; for the
  ///     closed set of metadata types it is `extern template` declared here and
  ///     instantiated in a dedicated translation unit.
  ///
  /// @return The built detray detector together with its name map.
  template <typename metadata_t>
  DetrayGeometry<metadata_t> convert(
      vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
      const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
      const std::string& detectorName = "") const;

 private:
  Config m_cfg;

  const Acts::Logger& logger() const { return *m_logger; }
  std::unique_ptr<const Acts::Logger> m_logger;
};

// Suppress implicit instantiation of `convert` for the closed set; the matching
// definitions are emitted in generated translation units (see CMakeLists.txt).
#define ACTS_DETRAY_EXTERN_CONVERT(META)                           \
  extern template DetrayGeometryConverter::DetrayGeometry<META>    \
  DetrayGeometryConverter::convert<META>(                          \
      vecmem::memory_resource&, const Acts::GeometryContext&,      \
      const std::shared_ptr<const Acts::TrackingGeometry>&,        \
      const std::string&) const;
ACTS_DETRAY_METADATA_FOR_EACH(ACTS_DETRAY_EXTERN_CONVERT)
#undef ACTS_DETRAY_EXTERN_CONVERT

}  // namespace ActsPlugins
