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
/// {func}`Acts::MagneticFieldProvider::makeCache`. This function accepts an
/// instance of {class}`Acts::MagneticFieldContext`, see
/// [](#magnetic-field-context) for details.
///
/// The main lookup method of {class}`Acts::MagneticFieldProvider` is @ref
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
/// // in event context
/// auto fieldContext = getExperimentFieldContext();
/// const Acts::MagneticFieldProvider& fieldProvider = getFieldProvider();
/// auto cache = fieldProvider.makeCache(fieldContext);
///
/// auto lookupResult = fieldProvider.getField(Acts::Vector3{10, 10, 10},
/// cache); if(!lookupResult.ok()) {
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
  using Cache = Acts::AnyBase<sizeof(char) * 512>;

  /// Make an opaque cache for the magnetic field. Instructs the specific
  /// implementation to generate a @c Cache instance for magnetic field lookup.
  ///
  /// @param mctx The magnetic field context to generate cache for
  /// @return Cache The opaque cache object
  virtual Cache makeCache(const MagneticFieldContext& mctx) const = 0;

  /// Retrieve magnetic field value at a given location. Requires a cache object
  /// created through makeCache().
  ///
  /// @param [in] position global 3D position for the lookup
  /// @param [in,out] cache Field provider specific cache object
  ///
  /// @return magnetic field vector at given position
  virtual Result<Vector3> getField(const Vector3& position,
                                   Cache& cache) const = 0;

  virtual ~MagneticFieldProvider() = default;
};

}  // namespace Acts
