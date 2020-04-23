// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

namespace BField {

/// The Context to be handed around
struct ScalableBFieldContext {
  double scalor = 1.;
};

/// @ingroup MagneticField
///
/// @brief returns a given constant field value at every point
///
/// This class is based on the constant magnetic field class
/// but allows a event based context
class ScalableBField final {
 public:
  struct Cache {
    double scalor = 1.;

    /// @brief constructor with context
    Cache(const Acts::MagneticFieldContext& mcfg) {
      scalor = std::any_cast<const ScalableBFieldContext>(mcfg).scalor;
    }
  };

  /// @brief construct constant magnetic field from field vector
  ///
  /// @param [in] B magnetic field vector in global coordinate system
  explicit ScalableBField(Acts::Vector3D B) : m_BField(std::move(B)) {}

  /// @brief construct constant magnetic field from components
  ///
  /// @param [in] Bx magnetic field component in global x-direction
  /// @param [in] By magnetic field component in global y-direction
  /// @param [in] Bz magnetic field component in global z-direction
  ScalableBField(double Bx = 0., double By = 0., double Bz = 0.)
      : m_BField(Bx, By, Bz) {}

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global position
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Acts::Vector3D getField(const Acts::Vector3D& /*position*/) const {
    return m_BField;
  }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global position
  /// @param [in] cache Cache object (is ignored)
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Acts::Vector3D getField(const Acts::Vector3D& /*position*/,
                          Cache& cache) const {
    return m_BField * cache.scalor;
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global position
  /// @param [out] derivative gradient of magnetic field vector as (3x3)
  /// matrix
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Acts::Vector3D getFieldGradient(
      const Acts::Vector3D& /*position*/,
      Acts::ActsMatrixD<3, 3>& /*derivative*/) const {
    return m_BField;
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global position
  /// @param [out] derivative gradient of magnetic field vector as (3x3)
  /// matrix
  /// @param [in] cache Cache object (is ignored)
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Acts::Vector3D getFieldGradient(const Acts::Vector3D& /*position*/,
                                  Acts::ActsMatrixD<3, 3>& /*derivative*/,
                                  Cache& cache) const {
    return m_BField * cache.scalor;
  }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  /// @note The method will always return true for the constant B-Field
  bool isInside(const Acts::Vector3D& /*position*/) const { return true; }

  /// @brief update magnetic field vector from components
  ///
  /// @param [in] Bx magnetic field component in global x-direction
  /// @param [in] By magnetic field component in global y-direction
  /// @param [in] Bz magnetic field component in global z-direction
  void setField(double Bx, double By, double Bz) { m_BField << Bx, By, Bz; }

  /// @brief update magnetic field vector
  ///
  /// @param [in] B magnetic field vector in global coordinate system
  void setField(const Acts::Vector3D& B) { m_BField = B; }

 private:
  /// magnetic field vector
  Acts::Vector3D m_BField;
};

}  // namespace BField
}  // namespace FW
