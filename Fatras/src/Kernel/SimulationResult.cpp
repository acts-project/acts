// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Kernel/SimulationResult.hpp"

#include <ostream>

std::ostream& ActsFatras::operator<<(std::ostream& os,
                                     SimulationParticleStatus status) {
  switch (status) {
    case SimulationParticleStatus::eAlive:
      return (os << "alive");
    case SimulationParticleStatus::eInteracted:
      return (os << "killed-by-interaction");
    case SimulationParticleStatus::eDecayed:
      return (os << "killed-by-decay");
    default:
      return (os << "invalid");
  }
  return os;
}
