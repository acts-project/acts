// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <functional>
#include <memory>

namespace Acts {

/// A field provider paired with its own reusable lookup cache.
/// Construct once per magnetic field context and reuse for sequential lookups.
/// The provider is kept alive by shared ownership. The context is borrowed:
/// it and any data referenced by its cache must outlive the query and remain
/// unchanged while the query is in use.
/// Do not share a query between concurrent operations; construct independent
/// queries for them (and for a new magnetic field context).
class MagneticFieldQuery {
 public:
  /// @param provider Field provider used for both cache creation and lookup
  /// @param mctx Borrowed context for cache creation and propagation
  /// @throws std::invalid_argument if provider is null
  MagneticFieldQuery(std::shared_ptr<const MagneticFieldProvider> provider,
                     const MagneticFieldContext& mctx);

  // A borrowed context must not be a temporary.
  MagneticFieldQuery(std::shared_ptr<const MagneticFieldProvider> provider,
                     const MagneticFieldContext&& mctx) = delete;

  MagneticFieldQuery(const MagneticFieldQuery&) = delete;
  MagneticFieldQuery& operator=(const MagneticFieldQuery&) = delete;
  /// Move another query, preserving its borrowed context.
  /// @param other Query to move from
  MagneticFieldQuery(MagneticFieldQuery&& other) = default;

  /// Move another query's provider, context reference and cache into this
  /// query.
  /// @param other Query to move from
  /// @return This query
  MagneticFieldQuery& operator=(MagneticFieldQuery&& other) = default;

  /// Access the context used to create this query's cache.
  /// @return Borrowed magnetic field context
  const MagneticFieldContext& context() const;

  /// Query the field, updating the private provider-specific cache.
  /// @param position Global position of the lookup
  /// @return Magnetic field or the provider's lookup error
  Result<Vector3> getField(const Vector3& position);

 private:
  std::shared_ptr<const MagneticFieldProvider> m_provider;
  std::reference_wrapper<const MagneticFieldContext> m_context;
  MagneticFieldProvider::Cache m_cache;
};

}  // namespace Acts
