// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::edm {

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool seed<BASE>::operator==(const seed<T>& other) const {
  return ((top_index() == other.top_index()) &&
          (middle_index() == other.middle_index()) &&
          (bottom_index() == other.bottom_index()));
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE std::strong_ordering seed<BASE>::operator<=>(
    const seed<T>& other) const {
  if (top_index() != other.top_index()) {
    return (top_index() <=> other.top_index());
  } else if (middle_index() != other.middle_index()) {
    return (middle_index() <=> other.middle_index());
  } else {
    return (bottom_index() <=> other.bottom_index());
  }
}

}  // namespace traccc::edm
