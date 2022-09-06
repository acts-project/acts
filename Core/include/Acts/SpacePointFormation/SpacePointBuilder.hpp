// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/SpacePointUtility.hpp"

#include <boost/container/static_vector.hpp>

namespace Acts {

/// @class SpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// measurements on the pixel or strip detectors need further treatment. This
/// class takes the measurements and provides the corresponding space points.
///
template <typename spacepoint_t>
class SpacePointBuilder {
 public:
  using Measurement = Acts::BoundVariantMeasurement;
  // Constructor
  /// @param cfg The configuration for the space point builder
  /// @param func The function that provides user's SP constructor with global pos, global cov, and sourceLinks.
  /// @param logger The logging instance
  SpacePointBuilder(SpacePointBuilderConfig cfg,
                    std::function<spacepoint_t(
                        Acts::Vector3, Acts::Vector2,
                        boost::container::static_vector<const SourceLink*, 2>)>
                        func,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("SpamcePointBuilder", Logging::INFO));

  // Default constructor
  SpacePointBuilder() = default;

  /// @brief Calculates the space points out of a given collection of measurements
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measurements vector of measurements
  /// @param opt option for the space point bulding. It contains the ends of the strips for strip SP building
  /// @param spacePointIt Output iterator for the space point
  template <template <typename...> typename container_t>
  void buildSpacePoint(
      const GeometryContext& gctx,
      const std::vector<const Measurement*>& measurements,
      const SpacePointBuilderOptions& opt,
      std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const;

  /// @brief Searches possible combinations of two measurements on different
  /// surfaces that may come from the same particles
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measurementsFront vector of measurements on a surface
  /// @param measurementsBack vector of measurements on another surface
  /// @param measurementPairs storage of the measurement pairs
  void makeMeasurementPairs(
      const GeometryContext& gctx,
      const std::vector<const Measurement*>& measurementsFront,
      const std::vector<const Measurement*>& measurementsBack,
      std::vector<std::pair<const Measurement*, const Measurement*>>&
          measurementPairs) const;

 protected:
  // configuration of the single hit space point builder
  SpacePointBuilderConfig m_config;

  /// @brief Function to create external space point
  /// The constructor of spacepoint_t with Vector3 global pos, Vector2 global
  /// cov, and vector of source link pointers.
  std::function<spacepoint_t(
      Acts::Vector3, Acts::Vector2,
      boost::container::static_vector<const SourceLink*, 2>)>
      m_spConstructor;

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  std::shared_ptr<const SpacePointUtility> m_spUtility;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
#include "Acts/SpacePointFormation/detail/SpacePointBuilder.ipp"
