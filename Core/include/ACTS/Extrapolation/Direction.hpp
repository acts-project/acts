// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIRECTION_HPP
#define ACTS_DIRECTION_HPP 1

namespace Acts {

/// @brief propagation direction relative to momentum
enum Direction : int { backward = -1, forward = 1 };

}  // namespace Acts
#endif  // ACTS_DIRECTION_HPP
