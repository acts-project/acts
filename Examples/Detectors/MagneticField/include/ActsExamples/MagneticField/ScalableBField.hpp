// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/BFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace ActsExamples {

/// The ScalableBField-specific magnetic field context.
struct ScalableBFieldContext {
  Acts::ActsScalar scalor = 1.;
};

/// A constant magnetic field that is scaled depending on the event context.
class ScalableBField final : public Acts::BFieldProvider {
 public:
  struct Cache {
    Acts::ActsScalar scalor = 1.;

    /// @brief constructor with context
    Cache(const Acts::MagneticFieldContext& mctx) {
      scalor = mctx.get<const ScalableBFieldContext>().scalor;
    }
  };

  /// @brief construct constant magnetic field from field vector
  ///
  /// @param [in] B magnetic field vector in global coordinate system
  explicit ScalableBField(Acts::Vector3 B) : m_BField(std::move(B)) {}

  /// @brief construct constant magnetic field from components
  ///
  /// @param [in] Bx magnetic field component in global x-direction
  /// @param [in] By magnetic field component in global y-direction
  /// @param [in] Bz magnetic field component in global z-direction
  ScalableBField(Acts::ActsScalar Bx = 0, Acts::ActsScalar By = 0,
                 Acts::ActsScalar Bz = 0)
      : m_BField(Bx, By, Bz) {}

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global position
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Acts::Vector3 getField(const Acts::Vector3& /*position*/) const override {
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
  Acts::Vector3 getField(const Acts::Vector3& /*position*/,
                         BFieldProvider::Cache& gCache) const override {
    Cache& cache = gCache.get<Cache>();
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
  Acts::Vector3 getFieldGradient(
      const Acts::Vector3& /*position*/,
      Acts::ActsMatrix<3, 3>& /*derivative*/) const override {
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
  Acts::Vector3 getFieldGradient(const Acts::Vector3& /*position*/,
                                 Acts::ActsMatrix<3, 3>& /*derivative*/,
                                 BFieldProvider::Cache& gCache) const override {
    Cache& cache = gCache.get<Cache>();
    return m_BField * cache.scalor;
  }

  Acts::BFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return Acts::BFieldProvider::Cache::make<Cache>(mctx);
  }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  /// @note The method will always return true for the constant B-Field
  bool isInside(const Acts::Vector3& /*position*/) const { return true; }

  /// @brief update magnetic field vector from components
  ///
  /// @param [in] Bx magnetic field component in global x-direction
  /// @param [in] By magnetic field component in global y-direction
  /// @param [in] Bz magnetic field component in global z-direction
  void setField(double Bx, double By, double Bz) { m_BField << Bx, By, Bz; }

  /// @brief update magnetic field vector
  ///
  /// @param [in] B magnetic field vector in global coordinate system
  void setField(const Acts::Vector3& B) { m_BField = B; }

 private:
  /// magnetic field vector
  Acts::Vector3 m_BField;
};  // namespace BField

}  // namespace ActsExamples
