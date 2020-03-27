// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <vector>

#include "Acts/EventData/MeasurementHelpers.hpp"

namespace FW {

/// A hit identifier with additional truth information.
///
/// In an addition to a unique identifier, an list of additional indices is
/// stored. These can refer e.g. to particle indices or truth hit indices.
/// Using indices instead of pointers allows more flexibility, i.e. we are not
/// fixed to a specific object type when using e.g. pointers, and is more
/// robust since no requirements on stable memory locations of the pointed-to
/// objects are necessary (as would be the case for pointers).
class SimIdentifier : public Acts::MinimalSourceLink {
 public:
  using Value = uint64_t;
  using Difference = int64_t;

  /// Constructor from encoded identifier value.
  ///
  /// @param value is the identifier value
  explicit SimIdentifier(Value value) : m_value(value) {}
  /// Constructor from encoded identifier value and truth information.
  ///
  /// @param value is the identifier value
  /// @param indices
  SimIdentifier(Value value, std::vector<std::size_t> indices)
      : m_value(value), m_indices(std::move(indices)) {}

  // Explicitely defaulted constructors and assignment operators
  SimIdentifier() = default;
  SimIdentifier(const SimIdentifier&) = default;
  SimIdentifier(SimIdentifier&&) = default;
  SimIdentifier& operator=(const SimIdentifier&) = default;
  SimIdentifier& operator=(SimIdentifier&&) = default;

  /// Assign from an identifier value.
  SimIdentifier& operator=(Value value);
  /// Cast to an identifier value.
  operator Value() const { return m_value; }

  /// Explicit access the underlying identifier value.
  Value value() const { return m_value; }
  /// Access all associated truth indices.
  const std::vector<std::size_t>& indices() const { return m_indices; }

  /// Attach a truth index to the identifier.
  void addIndex(std::size_t index) { m_indices.push_back(index); }

 private:
  /// The stored identifier value.
  Value m_value = 0u;
  /// Associated truth indices.
  std::vector<std::size_t> m_indices;

  friend constexpr bool operator<(const SimIdentifier& lhs,
                                  const SimIdentifier& rhs) {
    return lhs.m_value < rhs.m_value;
  }
  friend bool operator==(const SimIdentifier& lhs, const SimIdentifier& rhs) {
    return lhs.m_value == rhs.m_value;
  }
};

}  // end of namespace FW

using identifier_type = ::FW::SimIdentifier::Value;
using identifier_diff = ::FW::SimIdentifier::Difference;
using Identifier = ::FW::SimIdentifier;
