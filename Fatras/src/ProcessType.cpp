// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/ProcessType.hpp"

#include <ostream>

namespace ActsFatras {

std::ostream &operator<<(std::ostream &os, ProcessType processType) {
  os << static_cast<uint32_t>(processType);
  return os;
}

}  // namespace ActsFatras
