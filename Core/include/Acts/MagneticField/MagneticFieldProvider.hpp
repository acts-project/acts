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
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// Base class for all magnetic field providers
/// @ingroup magnetic_field
///
/// All magnetic field implementations inherit and implement from this
/// interface.
///
/// It provides a generic interface over different implementations. To speed up
/// magnetic field lookup, each implementation can have a *cache* object. The
/// cache object can for example be used to store a local interpolation domain,
/// to speed up nearby field lookups. The client is expected to pass into
/// lookup calls an instance of @ref Acts::MagneticFieldProvider::Cache.
///
/// The implementation is then free to use and update this cache instance as
/// needed. Before a client can issue field lookup calls, it needs to obtain an
/// initialized instance of this cache object. This can be achieved generically
/// for all implementations by using
/// @ref Acts::MagneticFieldProvider::makeCache. This function accepts an
/// instance of @ref Acts::MagneticFieldContext.
///
/// The main lookup method of @ref Acts::MagneticFieldProvider is @ref
/// Acts::MagneticFieldProvider::getField
///
/// Aside from the lookup position as a global position vector, it accepts an
/// instance of the opaque cache object mentioned before. The return value is a
/// @ref Acts::Result object. It either contains the field value at the
/// requested location, or an @ref Acts::MagneticFieldError in case of a lookup
/// failure, like an out-of-bounds lookup position.
///
/// Below is an example of how a client can interact with an instance of
/// @ref Acts::MagneticFieldProvider.
///
/// ```cpp
/// // In event context
/// auto fieldContext = getExperimentFieldContext();
/// const Acts::MagneticFieldProvider& fieldProvider = getFieldProvider();
/// // Make an opaque cache for field lookups
/// auto cache = fieldProvider.makeCache(fieldContext);
///
/// auto lookupResult = fieldProvider.getField(Acts::Vector3{10, 10, 10},
///                                   cache);
/// if(!lookupResult.ok()) {
///    throw std::runtime_error{"Field lookup failure"};
/// }
///
/// Acts::Vector3 fieldValue = *lookupResult;
/// ```
class MagneticFieldProvider {
 public:
  /// Opaque cache type that can store arbitrary implementation specific cache
  /// data. Examples are an interpolation cell, or an experiment specific
  /// conditions data handle.
  ///
  /// The cache is always creaded through @ref makeCache.
  using Cache = Acts::AnyBase<sizeof(char) * 512>;

  /// Magnetic field and its spatial gradient at one position
  struct FieldAndGradient {
    /// Magnetic field vector
    Vector3 field = Vector3::Zero();
    /// Spatial gradient of the field, with gradient(i, j) = dB_i / dx_j
    SquareMatrix3 gradient = SquareMatrix3::Zero();
  };

  /// Make an opaque cache for the magnetic field. Instructs the specific
  /// implementation to generate a @ref Acts::MagneticFieldProvider::Cache instance
  /// for magnetic field lookup.
  ///
  /// @param mctx The magnetic field context to generate cache for
  /// @return Cache The opaque cache object
  virtual Cache makeCache(const MagneticFieldContext& mctx) const = 0;

  /// Retrieve magnetic field value at a given location. Requires an instance
  /// of @ref Acts::MagneticFieldProvider::Cache created through @ref makeCache.
  ///
  /// @param [in] position global 3D position for the lookup
  /// @param [in,out] cache Field provider specific cache object
  ///
  /// @return magnetic field vector at given position
  virtual Result<Vector3> getField(const Vector3& position,
                                   Cache& cache) const = 0;

  /// Check if @ref getFieldAndGradient is implemented by this provider.
  ///
  /// A caller checks this once. If it is false, the caller can use
  /// @ref Acts::getFieldAndGradientNumerically instead.
  ///
  /// @return true if the provider implements the field gradient
  virtual bool providesFieldGradient() const { return false; }

  /// Retrieve the magnetic field and its spatial gradient at a given
  /// location.
  ///
  /// @param [in] position global 3D position for the lookup
  /// @param [in,out] cache Field provider specific cache object
  ///
  /// @return field and gradient at the given position, or
  ///         @ref MagneticFieldError::NotImplemented if
  ///         @ref providesFieldGradient is false
  virtual Result<FieldAndGradient> getFieldAndGradient(const Vector3& position,
                                                       Cache& cache) const;

  virtual ~MagneticFieldProvider() = default;
};

/// Compute the magnetic field gradient with central finite differences.
///
/// The function looks up the field at @p position and at the six positions
/// that are @p epsilon away along each axis. The error of the gradient is
/// proportional to @p epsilon squared for a smooth field. An interpolated
/// field map is only continuous at its cell edges, so there the result is an
/// average of the gradients of the neighbouring cells.
///
/// @param [in] provider The magnetic field provider
/// @param [in] position global 3D position for the lookup
/// @param [in,out] cache Field provider specific cache object
/// @param [in] epsilon Distance of the lookups from @p position
///
/// @return field and gradient at the given position, or the error of the
///         first lookup that failed
/// @ingroup magnetic_field
Result<MagneticFieldProvider::FieldAndGradient> getFieldAndGradientNumerically(
    const MagneticFieldProvider& provider, const Vector3& position,
    MagneticFieldProvider::Cache& cache, double epsilon);

}  // namespace Acts
