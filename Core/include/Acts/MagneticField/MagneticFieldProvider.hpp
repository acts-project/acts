// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/detail/SmallObjectCache.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <memory>

namespace Acts {

/// Base class for all magnetic field providers
class MagneticFieldProvider {
 public:
  using Cache = detail::SmallObjectCache;

  /// @brief Make an opaque cache for the magnetic field
  ///
  /// @param mctx The magnetic field context to generate cache for
  /// @return Cache The opaque cache object
  virtual Cache makeCache(const MagneticFieldContext& mctx) const = 0;

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  ///
  /// @return magnetic field vector at given position
  virtual Result<Vector3> getField(const Vector3& position,
                                   Cache& cache) const = 0;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @param [in,out] cache Cache object. Contains field cell used for
  /// @return magnetic field vector
  virtual Result<Vector3> getFieldGradient(const Vector3& position,
                                           ActsMatrix<3, 3>& derivative,
                                           Cache& cache) const = 0;

  virtual ~MagneticFieldProvider();
};

inline MagneticFieldProvider::~MagneticFieldProvider() = default;

}  // namespace Acts