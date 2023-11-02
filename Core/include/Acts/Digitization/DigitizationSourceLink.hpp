// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <vector>

namespace Acts {

/// Source link to connect digitization clusters back to truth information.
class DigitizationSourceLink final {
 public:
  /// Constructor from geometry identifier and truth indices.
  ///
  /// @param gid is the geometry identifier
  /// @param indices are the truth indices
  DigitizationSourceLink(GeometryIdentifier gid,
                         std::vector<std::size_t> indices = {})
      : m_geometryId(gid), m_indices(std::move(indices)) {}

  /// Construct and invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  DigitizationSourceLink() = default;
  DigitizationSourceLink(const DigitizationSourceLink&) = default;
  DigitizationSourceLink(DigitizationSourceLink&&) = default;
  DigitizationSourceLink& operator=(const DigitizationSourceLink&) = default;
  DigitizationSourceLink& operator=(DigitizationSourceLink&&) = default;

  /// Access all associated truth indices.
  const std::vector<std::size_t>& indices() const { return m_indices; }

  GeometryIdentifier geometryId() const { return m_geometryId; }

 private:
  GeometryIdentifier m_geometryId;
  /// Associated truth indices.
  std::vector<std::size_t> m_indices;
};

inline bool operator==(const DigitizationSourceLink& lhs,
                       const DigitizationSourceLink& rhs) {
  return (lhs.geometryId() == rhs.geometryId()) &&
         (lhs.indices() == rhs.indices());
}
inline bool operator!=(const DigitizationSourceLink& lhs,
                       const DigitizationSourceLink& rhs) {
  return !(lhs == rhs);
}

}  // namespace Acts
