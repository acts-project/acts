// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/RangeXD.hpp"

namespace Acts {

/// Magnetic field provider modelling a magnetic field consisting of
/// several (potentially overlapping) regions of constant values.
///
/// @ingroup magnetic_field
///
/// The multi-range constant field allows modelling cases where a magnetic field
/// can be described as multiple (potentially overlapping) regions, each of
/// which has its own constant magnetic field. This provides more flexibility
/// than the
/// @ref Acts::ConstantBField while providing higher performance than
/// @ref Acts::InterpolatedBFieldMap.
///
/// This magnetic field provider is configured using a list of pairs, where each
/// pair defines a region in three-dimensional space as well as a field vector.
/// Magnetic field lookup then proceeds by finding the *last* region in the
/// user-provided list that contains the requested coordinate and returning the
/// corresponding field vector.
///
/// The implementation uses a simple caching mechanism to store the last matched
/// region, providing improved performance for consecutive lookups within the
/// same region. This is thread-safe when each thread uses its own cache
/// instance. The field configuration itself is immutable after construction.
class MultiRangeBField final : public MagneticFieldProvider {
 private:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*unused*/);

    std::optional<std::size_t> index = {};
  };

  using BFieldRange = std::pair<RangeXD<3, double>, Vector3>;

  // The different ranges and their corresponding field vectors. Note that
  // regions positioned _later_ in this vector take priority over earlier
  // regions.
  std::vector<BFieldRange> fieldRanges;

 public:
  /// @brief Construct a magnetic field from a vector of ranges.
  ///
  /// @param ranges Vector of magnetic field ranges to use
  /// @warning These ranges are listed in increasing order of precedence,
  /// i.e. ranges further along the vector have higher priority.
  explicit MultiRangeBField(const std::vector<BFieldRange>& ranges);

  /// Construct from a vector of magnetic field ranges (move version).
  /// @param ranges Vector of magnetic field ranges to use (moved)
  /// @warning These ranges are listed in increasing order of precedence,
  /// i.e. ranges further along the vector have higher priority.
  explicit MultiRangeBField(std::vector<BFieldRange>&& ranges);

  /// @brief Construct a cache object.
  /// @param mctx Magnetic field context for cache creation
  /// @return Cache object for magnetic field computations
  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override;

  /// @brief Request the value of the magnetic field at a given position.
  ///
  /// @param [in] position Global 3D position for the lookup.
  /// @param [in, out] cache Cache object.
  /// @returns A successful value containing a field vector if the given
  /// location is contained inside any of the regions, or a failure value
  /// otherwise.
  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override;
};

}  // namespace Acts
