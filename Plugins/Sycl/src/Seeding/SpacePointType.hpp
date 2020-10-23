// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::Sycl::detail {

/// Type of a space point
enum class SpacePointType : int {
  Bottom = 0,  //< The referenced type is a "bottom" spacepoint
  Middle = 1,  //< The referenced type is a "middle" spacepoint
  Top = 2      //< The referenced type is a "top" spacepoint
};

}  // namespace Acts::Sycl::detail
