// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsPlugins/Root/TGeoAxes.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <tuple>

#include "RtypesCore.h"

class TGeoShape;
class TGeoMatrix;

namespace Acts {

class CylinderBounds;
class DiscBounds;
class PlanarBounds;
class Surface;

}  // namespace Acts

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

class TGeoDetectorElement;

/// Helper struct to convert TGeoShapes into Surface or Volume Bounds
struct TGeoSurfaceConverter {
  /// Convert a TGeoShape into cylinder surface components
  ///
  /// @param tgShape The TGeoShape
  /// @param rotation The rotation matrix as Double_t* from root
  /// @param translation The translation vector as Double_t* from root
  /// @param axes The axes definition
  /// @param scalor The unit scalor between TGeo and Acts
  ///
  /// @return tuple of DiscBounds, Transform, thickness
  static std::tuple<std::shared_ptr<const Acts::CylinderBounds>,
                    const Acts::Transform3, double>
  cylinderComponents(const TGeoShape& tgShape, const Double_t* rotation,
                     const Double_t* translation, TGeoAxes axes,
                     double scalor = 10.) noexcept(false);

  /// Convert a TGeoShape into disk surface components
  ///
  /// @param tgShape The TGeoShape
  /// @param rotation The rotation matrix as Double_t* from root
  /// @param translation The translation vector as Double_t* from root
  /// @param axes The axes definition
  /// @param scalor The unit scalor between TGeo and Acts
  ///
  /// @return tuple of DiscBounds, Transform, thickness
  static std::tuple<std::shared_ptr<const Acts::DiscBounds>,
                    const Acts::Transform3, double>
  discComponents(const TGeoShape& tgShape, const Double_t* rotation,
                 const Double_t* translation, TGeoAxes axes,
                 double scalor = 10.) noexcept(false);

  /// Convert a TGeoShape into plane surface components
  ///
  /// @param tgShape The TGeoShape
  /// @param rotation The rotation matrix as Double_t* from root
  /// @param translation The translation as a Double_t*
  /// @param axes The axes definition
  /// @param scalor The unit scalor between TGeo and Acts
  ///
  /// @return tuple of PlanarBounds, Transform, thickness
  static std::tuple<std::shared_ptr<const Acts::PlanarBounds>,
                    const Acts::Transform3, double>
  planeComponents(const TGeoShape& tgShape, const Double_t* rotation,
                  const Double_t* translation, TGeoAxes axes,
                  double scalor = 10.) noexcept(false);

  /// Convert a TGeoShape to a Surface
  ///
  /// @param tgShape The TGeoShape
  /// @param tgMatrix The matrix representing the transform
  /// @param axes The axes definition
  /// @param scalor The unit scalor between TGeo and Acts
  ///
  /// @return shared pointer to a surface and the original thickness that
  /// has been condensed to the surface
  static std::tuple<std::shared_ptr<Acts::Surface>, double> toSurface(
      const TGeoShape& tgShape, const TGeoMatrix& tgMatrix, TGeoAxes axes,
      double scalor = 10.) noexcept(false);

  /// Translate TGeo degree [0, 360) to radian
  /// * will correct to [-pi,pi)
  /// * it will return any multiple of 360.0 to 2pi
  /// @param deg The input in degree
  /// @return angle in radians
  static double toRadian(double deg) {
    double d = Acts::detail::wrap_periodic(deg, -180.0, 360.0);

    // If the angle was mapped to 0 from a non-0 multiple of 360, set it to 360
    if (constexpr double eps = 1e-6; std::abs(d) < eps && std::abs(deg) > eps) {
      d = 360.;
    }

    // Convert to rads
    return d * Acts::UnitConstants::degree;
  }
};

/// @}
}  // namespace ActsPlugins
