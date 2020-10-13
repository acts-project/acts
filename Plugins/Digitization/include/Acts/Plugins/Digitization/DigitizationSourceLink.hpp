// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cassert>
#include <vector>

namespace Acts {

class Surface;

/// Source link to connect digitization clusters back to truth information.
class DigitizationSourceLink final {
 public:
  /// Constructor from geometry identifier and truth indices.
  ///
  /// @param surface is the representing surface
  /// @param indices are the truth indices
  DigitizationSourceLink(const Surface& surface,
                         std::vector<std::size_t> indices = {})
      : m_surface(&surface), m_indices(std::move(indices)) {}

  /// Construct and invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  DigitizationSourceLink() = default;
  DigitizationSourceLink(const DigitizationSourceLink&) = default;
  DigitizationSourceLink(DigitizationSourceLink&&) = default;
  DigitizationSourceLink& operator=(const DigitizationSourceLink&) = default;
  DigitizationSourceLink& operator=(DigitizationSourceLink&&) = default;

  /// Access the reference surface.
  const Surface& referenceSurface() const {
    assert(m_surface and "Invalid Surface pointer in DigitizationSourceLink");
    return *m_surface;
  }
  /// Access all associated truth indices.
  const std::vector<std::size_t>& indices() const { return m_indices; }

 private:
  /// The stored surface. Use pointer to make the object copyable.
  const Surface* m_surface = nullptr;
  /// Associated truth indices.
  std::vector<std::size_t> m_indices;
};

inline bool operator==(const DigitizationSourceLink& lhs,
                       const DigitizationSourceLink& rhs) {
  return (lhs.referenceSurface() == rhs.referenceSurface()) and
         (lhs.indices() == rhs.indices());
}
inline bool operator!=(const DigitizationSourceLink& lhs,
                       const DigitizationSourceLink& rhs) {
  return not(lhs == rhs);
}

}  // namespace Acts
