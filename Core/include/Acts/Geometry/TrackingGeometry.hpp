// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

//  acts/Core/include/Acts/Geometry/TrackingGeometry.hpp
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/SurfaceVisitorConcept.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace Acts {

class Layer;
class Surface;
class PerigeeSurface;
class IMaterialDecorator;
class TrackingVolume;

using TrackingVolumePtr = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;

///  @class TrackingGeometry
///
///  The TrackingGeometry class is the owner of the constructed TrackingVolumes.
///
///  It enables both, a global search for an asociatedVolume
///  (respectively, if existing, a global search of an associated Layer or the
///  next associated Layer), such as a continuous navigation by BoundarySurfaces
///  between the confined TrackingVolumes.
class TrackingGeometry {
  /// Give the GeometryBuilder friend rights
  friend class TrackingGeometryBuilder;

 public:
  /// Constructor
  ///
  /// @param highestVolume is the world volume
  /// @param materialDecorator is a dediated decorator that can assign
  ///        surface or volume based material to the TrackingVolume
  /// @param hook Identifier hook to be applied to surfaces
  /// @param logger instance of a logger (defaulting to the "silent" one)
  TrackingGeometry(const MutableTrackingVolumePtr& highestVolume,
                   const IMaterialDecorator* materialDecorator = nullptr,
                   const GeometryIdentifierHook& hook = {},
                   const Logger& logger = getDummyLogger());

  /// Destructor
  ~TrackingGeometry();

  /// Access to the world volume
  /// @return plain pointer to the world volume
  const TrackingVolume* highestTrackingVolume() const;

  /// Access to the world volume
  /// @return shared pointer to the world volume
  const std::shared_ptr<const TrackingVolume>& highestTrackingVolumeShared()
      const;

  /// return the lowest tracking Volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to the lowest TrackingVolume
  const TrackingVolume* lowestTrackingVolume(const GeometryContext& gctx,
                                             const Vector3& gp) const;

  /// Forward the associated Layer information
  ///
  /// @param gctx is the context for this request (e.g. alignment)
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to assocaiated layer
  const Layer* associatedLayer(const GeometryContext& gctx,
                               const Vector3& gp) const;

  /// Register the beam tube
  ///
  /// @param beam is the beam line surface
  void registerBeamTube(std::shared_ptr<const PerigeeSurface> beam);

  /// @brief surface representing the beam pipe
  ///
  /// @note The ownership is not passed, e.g. do not delete the pointer
  ///
  /// @return raw pointer to surface representing the beam pipe
  ///         (could be a null pointer)
  const Surface* getBeamline() const;

  /// @brief Visit all sensitive surfaces
  ///
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  template <ACTS_CONCEPT(SurfaceVisitor) visitor_t>
  void visitSurfaces(visitor_t&& visitor) const {
    highestTrackingVolume()->template visitSurfaces<visitor_t>(
        std::forward<visitor_t>(visitor));
  }

  /// Search for a volume with the given identifier.
  ///
  /// @param id is the geometry identifier of the volume
  /// @retval nullptr if no such volume exists
  /// @retval pointer to the found volume otherwise.
  const TrackingVolume* findVolume(GeometryIdentifier id) const;

  /// Search for a surface with the given identifier.
  ///
  /// @param id is the geometry identifier of the surface
  /// @retval nullptr if no such surface exists
  /// @retval pointer to the found surface otherwise.
  const Surface* findSurface(GeometryIdentifier id) const;

  // adding new functions - implementation of misalignment of individual sensors

  /// Set the misalignment for a specific sensor (identified by  ID)  ///
  /// @param sensorId The identifier of the sensor
  /// @param misalignmentX misalignment in the x - direction
  /// @param misalignmentY misalignment in the y -  direction
  void setSensorMisalignment(GeometryIdentifier sensorId,
                             double misalignmentX, double misalignmentY);

  /// Get the misalignment value for a specific sensor identified by its ID
  ///
  /// @param sensorId The identifier of the sensor
  /// @return std::pair<double, double> The misalignment values in the
  ///         X and Y directions for the sensor
  std::pair<double, double> getSensorMisalignment(GeometryIdentifier sensorId) const;


  // Add functions to set and retrieve misalignment information
  void setSensorMisalignment(const GeometryIdentifier&, double misalignmentX, double misalignmentY);
  std::pair<double, double> getSensorMisalignment(const GeometryIdentifier& sensorName) const;

  void setCorrelatedMisalignment(double correlatedMisalignmentX, double correlatedMisalignmentY);
  std::pair<double, double> getCorrelatedMisalignment() const;


 private:
  // the known world
  TrackingVolumePtr m_world;
  // beam line
  std::shared_ptr<const PerigeeSurface> m_beam;
  // lookup containers
  std::unordered_map<GeometryIdentifier, const TrackingVolume*> m_volumesById;
  std::unordered_map<GeometryIdentifier, const Surface*> m_surfacesById;
  // misalignment map
  std::unordered_map<GeometryIdentifier, std::pair<double, double>> m_sensorMisalignment;
};

inline void TrackingGeometry::setSensorMisalignment(const GeometryIdentifier& sensor, double misalignmentX, double misalignmentY) {
  // Implementation to set the misalignment for the given sensor
  m_sensorMisalignment[sensor] = std::make_pair(misalignmentX, misalignmentY);
}


inline std::pair<double, double> TrackingGeometry::getSensorMisalignment(const GeometryIdentifier& sensor) const {
  // Implementation to retrieve the misalignment for the given sensor
  auto it = m_sensorMisalignment.find(sensor);
  if (it != m_sensorMisalignment.end()) {
    return it->second;
  } else {
    // Handle case when sensor misalignment is not found
    throw std::runtime_error("Misalignment information not available for sensor: "); // TODO
  }
}

inline void TrackingGeometry::setCorrelatedMisalignment(double correlatedMisalignmentX, double correlatedMisalignmentY) {
  // Implementation to set the correlated misalignment for all sensors
  for (auto& sensorMisalignment : m_sensorMisalignment) {
    sensorMisalignment.second = std::make_pair(correlatedMisalignmentX, correlatedMisalignmentY);
  }
}


inline std::pair<double, double> TrackingGeometry::getCorrelatedMisalignment() const {
  // We assume that all sensors have the same correlated misalignment
  // Thus, we can retrieve the misalignment from any sensor

  if (!m_sensorMisalignment.empty()) {
    // Retrieve the misalignment from the first sensor
    const auto& firstSensorMisalignment = m_sensorMisalignment.begin()->second;
    return firstSensorMisalignment;
  } else {
    // Return default values if no sensors are present
    return std::make_pair(0.0, 0.0);
  }
}

}  // namespace Acts

