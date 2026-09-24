// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"

namespace ActsExamples {

/// Output interface for finalized material-mapping results.
///
/// MaterialMapping passes the same TrackingGeometryMaterial to each writer in
/// its Config::materialWriters list, allowing one mapping run to produce
/// several outputs (for example, a versioned JSON map and a ROOT map). Each
/// writer owns its output configuration and serialization; the mapping
/// algorithm remains independent of the file format.
///
/// Writers are called once after map finalization, currently during destruction
/// of MaterialMapping. They do not receive per-event material tracks.
class IMaterialWriter {
 public:
  /// Virtual Destructor
  virtual ~IMaterialWriter() = default;

  /// Persist the finalized material assignments to the configured output.
  ///
  /// @param detMaterial Finalized assignments shared with all configured writers
  virtual void writeMaterial(
      const Acts::TrackingGeometryMaterial& detMaterial) = 0;
};

}  // namespace ActsExamples
