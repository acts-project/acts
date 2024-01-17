// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"

#if defined(__cpp_concepts)
#include <concepts>

namespace Acts {

template <typename T>
concept SurfaceVisitor = requires(T v) {
  {v(std::declval<const Surface*>())};
};

}  // namespace Acts

#endif
