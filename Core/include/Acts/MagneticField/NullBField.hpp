// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/BFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace Acts {

/// @ingroup MagneticField
/// @brief Null bfield which returns 0 always
class NullBField final : public BFieldProvider {
 public:
  struct Cache {
    /// @brief constructor with context
    Cache(const MagneticFieldContext& /*mcfg*/) {}
  };

  /// @brief Default constructor
  NullBField() = default;

  /// @copydoc BFieldBase::getField(const Vector3&)
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Vector3 getField(const Vector3& /*position*/) const override {
    return m_BField;
  }

  /// @copydoc BFieldBase::getField(const Vector3&,BFieldBase::Cache&)
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Vector3 getField(const Vector3& /*position*/,
                   BFieldProvider::Cache& /*cache*/) const override {
    return m_BField;
  }

  /// @copydoc BFieldBase::getFieldGradient(const Vector3&,ActsMatrix<3,3>&)
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3 getFieldGradient(const Vector3& /*position*/,
                           ActsMatrix<3, 3>& /*derivative*/) const override {
    return m_BField;
  }

  /// @copydoc BFieldBase::getFieldGradient(const
  /// Vector3&,ActsMatrix<3,3>&,BFieldBase::Cache&)
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3 getFieldGradient(const Vector3& /*position*/,
                           ActsMatrix<3, 3>& /*derivative*/,
                           BFieldProvider::Cache& /*cache*/) const override {
    return m_BField;
  }

  /// @copydoc BFieldBase::makeCache(const MagneticFieldContext&)
  Acts::BFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return Acts::BFieldProvider::Cache::make<Cache>(mctx);
  }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  /// @note The method will always return true for the null B-Field
  bool isInside(const Vector3& /*position*/) const { return true; }

 private:
  /// magnetic field vector
  const Vector3 m_BField = Vector3::Zero();
};
}  // namespace Acts
