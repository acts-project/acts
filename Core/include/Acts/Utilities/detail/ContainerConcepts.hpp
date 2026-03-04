// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail {

template <typename T>
concept ContainerHasAt = requires(const T &t) { t.at(0); };

template <typename T>
concept ContainerHasArrayAccess = requires(const T &t) { t[0]; };

}  // namespace Acts::detail
