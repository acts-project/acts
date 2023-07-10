// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

// acts/Core/include/Acts/Geometry/detail/DefaultDetectorElementBase.hpp 
namespace Acts {

class Surface;
// class GeometryIdentifier;

/// @class DetectorElementBase
///
/// This is the default base class for all tracking detector elements
/// with read-out relevant information. It provides the minimal interface
/// for the Acts proxy mechanism for surfaces, i.e. surfaces in the
/// Tracking geometry representing actual detection devices
///
class DetectorElementBase {
 public:
  DetectorElementBase() = default;
  virtual ~DetectorElementBase() = default;

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3& transform(const GeometryContext& gctx) const = 0;

  /// Return surface representation - const return pattern
  virtual const Surface& surface() const = 0;

  /// Non-const return pattern
  virtual Surface& surface() = 0;

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  virtual double thickness() const = 0;

  void setSensorMisalignment(GeometryIdentifier sensorId, double misalignmentX, double misalignmentY);
  std::pair<double, double> getSensorMisalignment(GeometryIdentifier sensorId) const;

  void setCorrelatedSensorMisalignment(GeometryIdentifier sensorId, double correlatedMisalignmentX, double correlatedMisalignmentY);
  std::pair<double, double> getCorrelatedSensorMisalignment(GeometryIdentifier sensorId) const;

 private:
  std::unordered_map<GeometryIdentifier, std::pair<double, double>> m_sensorMisalignment;
  std::unordered_map<GeometryIdentifier, std::pair<double, double>> m_correlatedSensorMisalignment;
};

void DetectorElementBase::setSensorMisalignment(GeometryIdentifier sensorId, double misalignmentX, double misalignmentY) {
  m_sensorMisalignment[sensorId] = std::make_pair(misalignmentX, misalignmentY);
}

std::pair<double, double> DetectorElementBase::getSensorMisalignment(GeometryIdentifier sensorId) const {
  auto it = m_sensorMisalignment.find(sensorId);
  if (it != m_sensorMisalignment.end()) {
    return it->second;
  }
  return std::make_pair(0.0, 0.0);  // Default values if sensor not found
}

void DetectorElementBase::setCorrelatedSensorMisalignment(GeometryIdentifier sensorId, double correlatedMisalignmentX, double correlatedMisalignmentY) {
  m_correlatedSensorMisalignment[sensorId] = std::make_pair(correlatedMisalignmentX, correlatedMisalignmentY);
}

std::pair<double, double> DetectorElementBase::getCorrelatedSensorMisalignment(GeometryIdentifier sensorId) const {
  auto it = m_correlatedSensorMisalignment.find(sensorId);
  if (it != m_correlatedSensorMisalignment.end()) {
    return it->second;
  }
  return std::make_pair(0.0, 0.0);  // Default values if sensor not found
}

}  // end of namespace Acts
