// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail {

template <template <typename...> class, template <typename...> class>
struct is_same_template : std::false_type {};

template <template <typename...> class T>
struct is_same_template<T, T> : std::true_type {};

}  // namespace Acts::detail
