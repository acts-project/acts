// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"

namespace Acts {
/// @brief SensitiveSurface is a surface which just overrides
///        the isSensitive check. It can be used to feed external
///        constraints into the track fit as free surfaces, i.e.
///        default GeometryIdentifier, without the need to also
///        construct a detector element
template <SurfaceConcept Surface_t>
class SensitiveSurface : public Surface_t {
 public:
  /// @brief Declare free surfaces to be sensitive
  bool isSensitive() const final { return true; }
  /// @brief  Declare the base class to be friend such
  ///         that the user can call
  ///         Surface::makeShared<SensitiveSurface<PlaneSurface>>();
  friend class Surface;
  /// @brief implement the assignment operators
  /// Assignment operator
  ///
  /// @param other is the source surface for copying
  /// @return Reference to this surface for assignment chaining
  SensitiveSurface& operator=(const SensitiveSurface& other) {
    static_cast<Surface_t&>(*this) = other;
    return *this;
  }
  /// @brief implement the assignment operators
  /// Assignment operator
  ///
  /// @param other is the source surface for copying
  /// @return Reference to this surface for assignment chaining
  SensitiveSurface& operator=(const Surface_t& other) {
    static_cast<Surface_t&>(*this) = other;
    return *this;
  }

 protected:
  /// @brief Copy the constructor from the base class. Enforce
  ///        the requirement that the base class inherits from Surface
  /// @param args: Base-class constructor arguments
  template <typename... args_t>
  explicit SensitiveSurface(args_t... args)
    requires(std::is_base_of_v<Surface, Surface_t>)
      : Surface_t(std::forward<args_t>(args)...) {}
};
}  // namespace Acts
