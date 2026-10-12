// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::details {

/// Type of a space point
enum class spacepoint_type : int {
  bottom = 0,  //< The referenced type is a "bottom" spacepoint
  middle = 1,  //< The referenced type is a "middle" spacepoint
  top = 2      //< The referenced type is a "top" spacepoint
};

}  // namespace traccc::details
