// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Set the Geometry Context PLUGIN
#ifdef ACTS_CORE_GEOMETRYCONTEXT_PLUGIN
#include ACTS_CORE_GEOMETRYCONTEXT_PLUGIN
#else

#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/detail/ContextType.hpp"

#include <ostream>

namespace Acts {

/// @brief This is the central definition of the Acts
/// payload object regarding detector geometry status (e.g. alignment)
///
/// @ingroup context
///
/// It is propagated through the code to allow for event/thread
/// dependent geometry changes.
///
/// ## Construction
///
/// For typical use cases with alignment or conditions data:
/// @code
/// MyAlignmentData alignment = ...;
/// GeometryContext gctx{alignment};
/// @endcode
///
/// For testing or simple applications without alignment:
/// @code
/// auto gctx = GeometryContext::dangerouslyDefaultConstruct();
/// @endcode
///
/// @note The default constructor is deprecated. Use the factory method
///       dangerouslyDefaultConstruct() to make empty context creation explicit.
class GeometryContext : public ContextType {
 public:
  /// Static factory method for default construction
  /// @note Use this when you need a default context for testing or
  ///       simple applications without alignment/conditions data
  static GeometryContext dangerouslyDefaultConstruct() {
    ACTS_PUSH_IGNORE_DEPRECATED()
    return GeometryContext();
    ACTS_POP_IGNORE_DEPRECATED()
  }

  /// Default constructor
  /// @deprecated Use GeometryContext::dangerouslyDefaultConstruct() instead
  ///             to make empty context construction explicit
  [[deprecated("Use GeometryContext::dangerouslyDefaultConstruct() instead")]]
  GeometryContext() = default;

  /// Move construct from arbitrary type (inherited from ContextType)
  /// @tparam T The type of the value to construct from
  /// @param value The value to construct from
  template <typename T>
    requires(!std::is_same_v<std::decay_t<T>, GeometryContext> &&
             !std::is_base_of_v<ContextType, std::decay_t<T> >)
  explicit GeometryContext(T&& value) : ContextType(std::forward<T>(value)) {}

  /// Copy construct from arbitrary type (inherited from ContextType)
  /// @tparam T The type of the value to construct from
  /// @param value The value to construct from
  template <typename T>
    requires(!std::is_same_v<std::decay_t<T>, GeometryContext> &&
             !std::is_base_of_v<ContextType, std::decay_t<T> >)
  explicit GeometryContext(const T& value) : ContextType(value) {}

  /// Inherit assignment operators
  using ContextType::operator=;
};

/// Helper struct that stores an object and a context, and will print it to
/// an outstream.
/// This allows you to write
/// ```cpp
/// std::cout << surface->toStream(geoContext) << std::endl;
/// ```
template <typename T>
struct GeometryContextOstreamWrapper {
  /// Constructor wrapping an object with geometry context for output streaming
  /// @param object The object to wrap for streaming
  /// @param gctx The geometry context to associate with the object
  GeometryContextOstreamWrapper(const T& object, const GeometryContext& gctx)
      : m_object(object), m_gctx(gctx) {}

  GeometryContextOstreamWrapper(const GeometryContextOstreamWrapper&) = delete;
  GeometryContextOstreamWrapper(GeometryContextOstreamWrapper&&) = delete;
  GeometryContextOstreamWrapper& operator=(
      const GeometryContextOstreamWrapper&) = delete;
  GeometryContextOstreamWrapper& operator=(GeometryContextOstreamWrapper&&) =
      delete;

  friend std::ostream& operator<<(
      std::ostream& os, const GeometryContextOstreamWrapper& osWrapper) {
    osWrapper.toStream(os);
    return os;
  }

 private:
  void toStream(std::ostream& os) const { m_object.toStreamImpl(m_gctx, os); }

  const T& m_object;
  const GeometryContext& m_gctx;
};

}  // namespace Acts

#endif
