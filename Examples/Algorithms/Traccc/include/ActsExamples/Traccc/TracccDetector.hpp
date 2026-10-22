// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

#include <traccc/bfield/magnetic_field.hpp>
#include <traccc/geometry/detector.hpp>
#include <traccc/geometry/detector_buffer.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsExamples {

/// Detector and magnetic field set up once for a traccc device backend
///
/// Holds the (in-memory converted) detray detector of the ODD on the host
/// and its copy on the device, the host and device magnetic field, and the
/// mapping between detray and ACTS surface identifiers. It is shared between
/// all traccc device algorithms of a sequence.
class TracccDetector {
 public:
  /// The detray detector type used for the ODD
  using DetectorType = traccc::odd_detector;

  struct Config {
    /// The backend to set the detector up for
    std::shared_ptr<const TracccBackend> backend;
    /// The ACTS tracking geometry the detray detector was converted from
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Host memory resource the detray detector was created with. It has to
    /// outlive the detector, so it is declared (and destroyed) before it.
    std::shared_ptr<vecmem::memory_resource> detrayMemoryResource;
    /// The detray detector converted from the ACTS tracking geometry
    std::shared_ptr<const DetectorType> detrayDetector;
    /// Constant magnetic field, used if no field file is given
    Acts::Vector3 constantBField = Acts::Vector3::Zero();
    /// Optional covfie magnetic field file (overrides the constant field)
    std::string bFieldFile;
  };

  explicit TracccDetector(const Config& cfg,
                          std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  const TracccBackend& backend() const { return *m_cfg.backend; }
  const Acts::TrackingGeometry& trackingGeometry() const {
    return *m_cfg.trackingGeometry;
  }
  const DetectorType& hostDetector() const { return *m_cfg.detrayDetector; }
  const traccc::detector_buffer& deviceDetector() const {
    return m_deviceDetector;
  }

  const traccc::magnetic_field& hostField() const { return m_hostField; }
  const traccc::magnetic_field& deviceField() const { return m_deviceField; }

  /// ACTS geometry identifier of a detray surface (sensitive surfaces only)
  std::optional<Acts::GeometryIdentifier> actsIdentifier(
      std::uint64_t detrayId) const;
  /// detray surface identifier of an ACTS surface (sensitive surfaces only)
  std::optional<std::uint64_t> detrayIdentifier(
      Acts::GeometryIdentifier actsId) const;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  traccc::detector_buffer m_deviceDetector;
  traccc::magnetic_field m_hostField;
  traccc::magnetic_field m_deviceField;

  std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> m_detrayToActs;
  std::unordered_map<Acts::GeometryIdentifier, std::uint64_t> m_actsToDetray;
};

}  // namespace ActsExamples
