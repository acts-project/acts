// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @ingroup MagneticField
///
/// @brief allows to use a shared magnetic field
/// in several places and with multiple steppers
/// mainly targeted to save memory
template <typename BField>
class SharedBField {
 public:
  // typedef wrapped BField's cache type
  using Cache = typename BField::Cache;

  /// Constructur with a shared pointer from a shared pointer.
  /// @note Since it is a shared field, we enforce it to be const
  /// @tparam bField is the shared BField to be stored
  SharedBField(std::shared_ptr<const BField> bField) : m_bField(bField) {}

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  ///
  /// @return magnetic field vector at given position
  Vector3D getField(const Vector3D& position) const {
    return m_bField->getField(position);
  }

  /// @brief Retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  /// @param [in,out] cache Cache object, passed through to wrapped BField
  Vector3D getField(const Vector3D& position, Cache& cache) const {
    return m_bField->getField(position, cache);
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D getFieldGradient(const Vector3D& position,
                            ActsMatrixD<3, 3>& derivative) const {
    return m_bField->getFieldGradient(position, derivative);
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @param [in,out] cache Cache object, passed through to wrapped BField
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D getFieldGradient(const Vector3D& position,
                            ActsMatrixD<3, 3>& derivative, Cache& cache) const {
    return m_bField->getFieldGradient(position, derivative, cache);
  }

 private:
  std::shared_ptr<const BField> m_bField;
};

}  // namespace Acts
