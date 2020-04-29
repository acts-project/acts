// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include <vector>
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class ISurfaceMaterial
///
/// Virtual base class of surface based material description
///
/// MaterialProperties that are associated to a surface,
/// extended by certain special representations (binned, homogenous)
///
class ISurfaceMaterial {
 public:
  /// Constructor
  ISurfaceMaterial() = default;

  /// Constructor
  ///
  /// @param splitFactor is the splitting ratio between pre/post update
  ISurfaceMaterial(double splitFactor) : m_splitFactor(splitFactor) {}

  /// Destructor
  virtual ~ISurfaceMaterial() = default;

  /// Scale operator
  ///
  /// @param scale is the scale factor applied
  virtual ISurfaceMaterial& operator*=(double scale) = 0;

  /// Return method for full material description of the Surface
  /// - from local coordinate on the surface
  ///
  /// @param lp is the local position used for the (eventual) lookup
  ///
  /// @return const MaterialProperties
  virtual const MaterialProperties& materialProperties(
      const Vector2D& lp) const = 0;

  /// Return method for full material description of the Surface
  /// - from the global coordinates
  ///
  /// @param gp is the global position used for the (eventual) lookup
  ///
  /// @return const MaterialProperties
  virtual const MaterialProperties& materialProperties(
      const Vector3D& gp) const = 0;

  /// Direct access via bins to the MaterialProperties
  ///
  /// @param ib0 is the material bin in dimension 0
  /// @param ib1 is the material bin in dimension 1
  virtual const MaterialProperties& materialProperties(size_t ib0,
                                                       size_t ib1) const = 0;

  /// Update pre factor
  ///
  /// @param pDir is the navigation direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  double factor(NavigationDirection pDir, MaterialUpdateStage mStage) const;

  /// Return method for fully scaled material description of the Surface
  /// - from local coordinate on the surface
  ///
  /// @param lp is the local position used for the (eventual) lookup
  /// @param pDir is the navigation direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  ///
  /// @return MaterialProperties
  MaterialProperties materialProperties(const Vector2D& lp,
                                        NavigationDirection pDir,
                                        MaterialUpdateStage mStage) const;

  /// Return method for full material description of the Surface
  /// - from the global coordinates
  ///
  /// @param gp is the global position used for the (eventual) lookup
  /// @param pDir is the navigation direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  ///
  /// @return MaterialProperties
  MaterialProperties materialProperties(const Vector3D& gp,
                                        NavigationDirection pDir,
                                        MaterialUpdateStage mStage) const;

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual ISurfaceMaterial::toStream method
  ///
  /// @return modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const ISurfaceMaterial& sm) {
    sm.toStream(out);
    return out;
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  virtual std::ostream& toStream(std::ostream& sl) const = 0;

 protected:
  double m_splitFactor{1.};  //!< the split factor in favour of oppositePre
};

inline double ISurfaceMaterial::factor(NavigationDirection pDir,
                                       MaterialUpdateStage mStage) const {
  if (mStage == Acts::fullUpdate) {
    return 1.;
  }
  return (pDir * mStage > 0 ? m_splitFactor : 1. - m_splitFactor);
}

inline MaterialProperties ISurfaceMaterial::materialProperties(
    const Vector2D& lp, NavigationDirection pDir,
    MaterialUpdateStage mStage) const {
  // The plain material properties associated to this bin
  MaterialProperties plainMatProp = materialProperties(lp);
  // Scale if you have material to scale
  if (plainMatProp) {
    double scaleFactor = factor(pDir, mStage);
    if (scaleFactor == 0.) {
      return MaterialProperties();
    }
    plainMatProp.scaleThickness(scaleFactor);
  }
  return plainMatProp;
}

inline MaterialProperties ISurfaceMaterial::materialProperties(
    const Vector3D& gp, NavigationDirection pDir,
    MaterialUpdateStage mStage) const {
  // The plain material properties associated to this bin
  MaterialProperties plainMatProp = materialProperties(gp);
  // Scale if you have material to scale
  if (plainMatProp) {
    double scaleFactor = factor(pDir, mStage);
    if (scaleFactor == 0.) {
      return MaterialProperties();
    }
    plainMatProp.scaleThickness(scaleFactor);
  }
  return plainMatProp;
}

}  // namespace Acts
