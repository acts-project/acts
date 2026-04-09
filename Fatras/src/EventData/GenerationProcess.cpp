// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/GenerationProcess.hpp"

#include <ostream>

namespace ActsFatras {

std::ostream &operator<<(std::ostream &os, GenerationProcess processType) {
  using enum GenerationProcess;

  switch (processType) {
    case eUndefined:
      return (os << "undefined");
    case eDecay:
      return (os << "decay");
    case ePhotonConversion:
      return (os << "photon conversion");
    case eBremsstrahlung:
      return (os << "bremsstrahlung");
    case eNuclearInteraction:
      return (os << "nuclear interaction");
    default:
      return (os << static_cast<std::uint32_t>(processType));
  }
}

}  // namespace ActsFatras
