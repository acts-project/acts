// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsFatras/EventData/Hit.hpp"

namespace FW {

/// Source link class for simulation in the acts-framework.
///
/// The source link stores the measuremts, surface, and the associated simulated
/// truth hit.
///
/// @todo Allow multiple truth hits e.g. for merged hits.
class SimSourceLink {
 public:
  SimSourceLink(const Acts::Surface& surface, const ActsFatras::Hit& truthHit,
                size_t dim, Acts::BoundVector values, Acts::BoundMatrix cov)
      : m_values(values),
        m_cov(cov),
        m_dim(dim),
        m_geometryId(truthHit.geometryId()),
        m_surface(&surface),
        m_truthHit(&truthHit) {}
  /// Must be default_constructible to satisfy SourceLinkConcept.
  SimSourceLink() = default;
  SimSourceLink(SimSourceLink&&) = default;
  SimSourceLink(const SimSourceLink&) = default;
  SimSourceLink& operator=(SimSourceLink&&) = default;
  SimSourceLink& operator=(const SimSourceLink&) = default;

  constexpr Acts::GeometryID geometryId() const { return m_geometryId; }
  constexpr const Acts::Surface& referenceSurface() const { return *m_surface; }
  constexpr const ActsFatras::Hit& truthHit() const { return *m_truthHit; }

  Acts::FittableMeasurement<SimSourceLink> operator*() const {
    if (m_dim == 0) {
      throw std::runtime_error("Cannot create dim 0 measurement");
    } else if (m_dim == 1) {
      return Acts::Measurement<SimSourceLink, Acts::ParDef::eLOC_0>{
          m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<1, 1>(),
          m_values[0]};
    } else if (m_dim == 2) {
      return Acts::Measurement<SimSourceLink, Acts::ParDef::eLOC_0,
                               Acts::ParDef::eLOC_1>{
          m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
          m_values[0], m_values[1]};
    } else {
      throw std::runtime_error("Dim " + std::to_string(m_dim) +
                               " currently not supported.");
    }
  }

 private:
  Acts::BoundVector m_values;
  Acts::BoundMatrix m_cov;
  size_t m_dim = 0u;
  // store geo id copy to avoid indirection via truth hit
  Acts::GeometryID m_geometryId;
  // need to store pointers to make the object copyable
  const Acts::Surface* m_surface;
  const ActsFatras::Hit* m_truthHit;

  friend constexpr bool operator==(const SimSourceLink& lhs,
                                   const SimSourceLink& rhs) {
    return lhs.m_truthHit == rhs.m_truthHit;
  }
};

/// Store source links ordered by geometry identifier.
using SimSourceLinkContainer = GeometryIdMultiset<SimSourceLink>;

}  // end of namespace FW
