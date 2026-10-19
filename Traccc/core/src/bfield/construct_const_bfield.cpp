// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/bfield/construct_const_bfield.hpp"

namespace traccc {

magnetic_field construct_const_bfield(const vector3& v) {
  return construct_const_bfield(v[0], v[1], v[2]);
}

}  // namespace traccc
