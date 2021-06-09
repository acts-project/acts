// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <memory>

namespace Acts {

/// @ingroup MagneticField
///
/// @brief allows to use a shared magnetic field
/// in several places and with multiple steppers
/// mainly targeted to save memory
template <typename BField>
class SharedBField final : public MagneticFieldProvider {
 public:
  using Cache = typename BField::Cache;

  /// Disallow construction without a valid underlying field.
  SharedBField() = delete;

  /// Constructur with a shared pointer from a shared pointer.
  /// @note Since it is a shared field, we enforce it to be const
  /// @tparam bField is the shared BField to be stored
  SharedBField(std::shared_ptr<const BField> bField) : m_bField(bField) {}

  /// @copydoc MagneticFieldProvider::getField(const
  /// Vector3&,MagneticFieldProvider::Cache&)
  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override {
    return m_bField->getField(position, cache);
  }

  /// @copydoc MagneticFieldProvider::getFieldGradient(const
  /// Vector3&,ActsMatrix<3,3>&,MagneticFieldProvider::Cache&)
  Result<Vector3> getFieldGradient(
      const Vector3& position, ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override {
    return m_bField->getFieldGradient(position, derivative, cache);
  }

  /// @copydoc MagneticFieldProvider::makeCache(const MagneticFieldContext&)
  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return m_bField->makeCache(mctx);
  }

 private:
  std::shared_ptr<const BField> m_bField;
};

}  // namespace Acts
