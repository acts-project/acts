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
#include "Acts/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

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
  /// @param cfg the configuration for the space point builder
  /// @param logger The logging instance
  SpacePointBuilder(SpacePointBuilderConfig cfg,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("SpacePointBuilder", Logging::INFO));

  // Default constructor
  SpacePointBuilder() = default;

  /// @brief Calculates the space points out of a given collection of measurements
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param spacePointIt Output iterator for the space points
  /// @param frontMeasurements measurements on the front surfaces for the strip SP formation, or all measurements for the pixel SP.
  /// @param backMeasurements measurements on the back surfaces for the strip SP. nullptr for pixel SP formation
  template <template <typename...> typename container_t>
  void calculateSpacePoints(
      const GeometryContext& gctx,
      std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt,
      const std::vector<const Measurement*>* frontMeasurements,
      const std::vector<const Measurement*>* backMeasurements = nullptr) const;

 protected:
  /// @brief Getter method for the local coordinates of a measurement
  /// on its corresponding surface
  ///
  /// @param meas measurement that holds the neccesary information of the hit position.
  /// @return vector of the local coordinates of the measurement on the surface
  Vector2 getLocalPos(const Measurement& meas) const;
  std::pair<Acts::Vector2, Acts::SymMatrix2> getLocalPosCov(
      const Measurement& meas) const;

  /// @brief Getter method for the global coordinates of a measurement
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param meas measurement that holds the necessary
  /// information
  /// @return vectors of the global coordinates and covariance of the measurement
  std::pair<Vector3, Vector2> globalCoords(const GeometryContext& gctx,
                                           const Measurement& meas) const;

  /// @brief Calculates the space points out of a given collection of measurements
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measurements vector of measurements
  /// @param spacePointIt Output iterator for the space points
  template <template <typename...> typename container_t>
  void calculateSingleHitSpacePoints(
      const GeometryContext& gctx,
      const std::vector<const Measurement*>& measurements,
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

  /// @brief Searches possible combinations of two measurements on different
  /// surfaces that may come from the same particles
  /// @param gctx The geometry context to use
  /// @param measurementPairs pairs of measurements that are space point candidates
  /// @param spacePointIt storage of the results
  /// @note If no configuration is set, the default values will be used
  template <template <typename...> typename container_t>
  void calculateDoubleHitSpacePoints(
      const Acts::GeometryContext& gctx,
      const std::vector<std::pair<const Measurement*, const Measurement*>>&
          measurementPairs,
      std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const;

  /// @brief Calculates the top and bottom ends of a strip detector element
  /// that corresponds to a given hit
  /// @param gctx The geometry context to use
  /// @param measurement object that stores the information about the hit
  /// @return vectors to the top and bottom end of the SDE
  std::pair<Acts::Vector3, Acts::Vector3> endsOfStrip(
      const Acts::GeometryContext& gctx, const Measurement& measurement) const;

  /// @brief Get global covariance from the local position and covariance
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param geoId The geometry ID
  /// @param localPos The local position
  /// @param localCov The local covariance matrix
  /// @return (rho, z) components of the global covariance
  Acts::Vector2 globalCov(const Acts::GeometryContext& gctx,
                          const Acts::GeometryIdentifier& geoId,
                          const Acts::Vector2& localPos,
                          const Acts::SymMatrix2& localCov) const;

  /// @brief Get the first component of the local covariance.
  /// @param meas The measurement
  /// @return the (0, 0) component of the local covariance
  double getLoc0Var(const Measurement& meas) const;

  /// @brief Calculate the global covariance from the front and back measurement in the strip SP formation
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measFront The measurement on the front layer
  /// @param measBack The measurement on the back layer
  /// @param theta The angle between the two strips
  /// @return (rho, z) components of the global covariance
  Acts::Vector2 calcGlobalVars(const Acts::GeometryContext& gctx,
                               const Measurement& measFront,
                               const Measurement& measBack,
                               const double theta) const;

  /// @brief Get source link from the measurement
  /// @param meas The measurement
  const Acts::SourceLink* getSourceLink(const Measurement meas) const;

  // configuration of the single hit space point builder
  SpacePointBuilderConfig m_config;

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
#include "Acts/SpacePointFormation/detail/SpacePointBuilder.ipp"
