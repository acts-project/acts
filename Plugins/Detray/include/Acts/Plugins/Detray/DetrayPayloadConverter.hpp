// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <memory>

namespace detray::io {
struct detector_payload;
struct transform_payload;
struct mask_payload;
}  // namespace detray::io

namespace Acts {

class GeometryContext;
class TrackingGeometry;
class SurfaceBounds;

class DetrayPayloadConverter {
 public:
  static detray::io::detector_payload convert(
      const GeometryContext& gctx, const TrackingGeometry& trackingGeometry);

  static detray::io::transform_payload convertTransform(
      const Transform3& transform);

  /// @param forPortal detray special cases the local parametrization for portals for performance reasons
  static detray::io::mask_payload convertMask(const Acts::SurfaceBounds& bounds,
                                              bool forPortal);

 private:
  explicit DetrayPayloadConverter(const GeometryContext& gctx,
                                  const TrackingGeometry& trackingGeometry);

  ~DetrayPayloadConverter();

  void convertPayload();

  std::unique_ptr<detray::io::detector_payload> m_detectorPayload;
};
}  // namespace Acts
