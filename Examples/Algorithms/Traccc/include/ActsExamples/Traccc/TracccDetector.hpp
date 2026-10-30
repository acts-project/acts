// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"

#include <memory>

#include <traccc/geometry/detector_buffer.hpp>
#include <traccc/geometry/host_detector.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsExamples {

/// Detector set up once for a traccc device backend
///
/// Holds the (in-memory converted) detray detector on the host and its copy
/// on the device, and the mapping between detray and ACTS surface
/// identifiers that was created by the conversion. It is shared between all
/// traccc device algorithms of a sequence. The detray detector is type erased,
/// so any detector type that traccc was built with is supported.
class TracccDetector {
 public:
  struct Config {
    /// The backend to set the detector up for
    std::shared_ptr<const TracccBackend> backend;
    /// The ACTS tracking geometry the detray detector was converted from
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Host memory resource the detray detector was created with. It has to
    /// outlive the detector, so it is declared (and destroyed) before it.
    std::shared_ptr<vecmem::memory_resource> detrayMemoryResource;
    /// The detray detector converted from the ACTS tracking geometry
    std::shared_ptr<const traccc::host_detector> detrayDetector;
    /// The ACTS <-> detray surface identifier map from the same conversion
    std::shared_ptr<const ActsPlugins::DetraySurfaceIdMap> surfaceIds;
  };

  TracccDetector(const Config& cfg, std::unique_ptr<const Acts::Logger> logger);

  const Config& config() const { return m_cfg; }

  const TracccBackend& backend() const { return *m_cfg.backend; }
  const Acts::TrackingGeometry& trackingGeometry() const {
    return *m_cfg.trackingGeometry;
  }
  const traccc::host_detector& hostDetector() const {
    return *m_cfg.detrayDetector;
  }
  const traccc::detector_buffer& deviceDetector() const {
    return m_deviceDetector;
  }

  /// ACTS <-> detray identifier map of the sensitive surfaces
  const ActsPlugins::DetraySurfaceIdMap& surfaceIds() const {
    return *m_cfg.surfaceIds;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  traccc::detector_buffer m_deviceDetector;
};

}  // namespace ActsExamples
