// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/make_magnetic_field.hpp"

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/io/read_magnetic_field.hpp"

namespace traccc::details {

magnetic_field make_magnetic_field(const opts::magnetic_field& opts) {
  if (opts.read_from_file) {
    magnetic_field result;
    io::read_magnetic_field(result, opts.file, opts.format);
    return result;
  } else {
    return construct_const_bfield({0.f, 0.f, opts.value});
  }
}

}  // namespace traccc::details
