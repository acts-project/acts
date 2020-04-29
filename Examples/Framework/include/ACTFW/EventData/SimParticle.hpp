// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/container/flat_set.hpp>

#include "ActsFatras/EventData/Particle.hpp"

namespace FW {
namespace detail {
struct CompareParticleId {
  using is_transparent = void;
  constexpr bool operator()(const ActsFatras::Particle& lhs,
                            const ActsFatras::Particle& rhs) const {
    return lhs.particleId() < rhs.particleId();
  }
  constexpr bool operator()(ActsFatras::Barcode lhs,
                            const ActsFatras::Particle& rhs) const {
    return lhs < rhs.particleId();
  }
  constexpr bool operator()(const ActsFatras::Particle& lhs,
                            ActsFatras::Barcode rhs) const {
    return lhs.particleId() < rhs;
  }
};
}  // namespace detail

using SimParticle = ::ActsFatras::Particle;
/// Store particles ordered by particle identifier.
using SimParticleContainer =
    ::boost::container::flat_set<::ActsFatras::Particle,
                                 detail::CompareParticleId>;

}  // end of namespace FW
