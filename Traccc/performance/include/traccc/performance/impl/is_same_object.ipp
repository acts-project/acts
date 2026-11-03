// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::details {

template <typename T>
is_same_object<T>::is_same_object(const T& ref, scalar) : m_ref(ref) {}

template <typename T>
bool is_same_object<T>::operator()(const T& obj) const {
  return (obj == m_ref.get());
}

}  // namespace traccc::details
