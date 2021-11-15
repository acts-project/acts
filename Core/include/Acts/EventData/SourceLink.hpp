// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

namespace Acts {

/// Base class for all SourceLink objects. Exposes a minimal nonvirtual
/// interface
class SourceLink {
 protected:
  /// Constructor for @c SourceLink. Protected so it cannot be instantiated on it's own
  /// @param id The geometry identifier this source link is associated with
  constexpr SourceLink(GeometryIdentifier id) : m_geometryId{id} {}

 public:
  /// Getter for the geometry identifier
  /// @return The GeometryIdentifier
  constexpr GeometryIdentifier geometryId() const { return m_geometryId; }

  /// Virtual destructor, required for safely storing source links as their base
  virtual ~SourceLink() = 0;

 private:
  GeometryIdentifier m_geometryId;
};

inline SourceLink::~SourceLink() = default;

}  // namespace Acts
