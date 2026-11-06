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

namespace Acts {

/// A field provider paired with its own reusable lookup cache.
/// Construct once per magnetic field context and reuse for sequential lookups.
/// The provider and context are borrowed and must outlive the query. The
/// context and any context-dependent data referenced by the cache must remain
/// unchanged while the query is in use. Do not share a query between concurrent
/// operations; construct independent queries for them (and for a new magnetic
/// field context).
class MagneticFieldQuery {
 public:
  /// @param provider Field provider used for both cache creation and lookup
  /// @param mctx Borrowed context for cache creation and propagation
  MagneticFieldQuery(const MagneticFieldProvider& provider,
                     const MagneticFieldContext& mctx);

  // Borrowed objects must not be temporaries.
  MagneticFieldQuery(const MagneticFieldProvider&& provider,
                     const MagneticFieldContext& mctx) = delete;

  MagneticFieldQuery(const MagneticFieldProvider& provider,
                     const MagneticFieldContext&& mctx) = delete;

  /// Copy the cache, retaining the borrowed provider and context.
  /// @param other Query whose cache state is copied
  MagneticFieldQuery(const MagneticFieldQuery& other) = default;

  /// Copy the cache and borrowed provider/context references.
  /// @param other Query whose cache state is copied
  /// @return This query
  MagneticFieldQuery& operator=(const MagneticFieldQuery& other) = default;
  /// Move another query, preserving its borrowed provider and context.
  /// @param other Query to move from
  MagneticFieldQuery(MagneticFieldQuery&& other) = default;

  /// Move another query's provider and context references and cache into this
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
  // Transitional access for tools accepting legacy cache-taking linearizers.
  // Remove with their legacy configuration members.
  friend class FullBilloirVertexFitter;
  friend class AdaptiveMultiVertexFitter;
  friend class IterativeVertexFinder;

  std::reference_wrapper<const MagneticFieldProvider> m_provider;
  std::reference_wrapper<const MagneticFieldContext> m_context;
  MagneticFieldProvider::Cache m_cache;
};

}  // namespace Acts
