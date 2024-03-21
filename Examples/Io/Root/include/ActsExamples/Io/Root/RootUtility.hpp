// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>

#include <TMathBase.h>

namespace ActsExamples {
namespace RootUtility {

/// Sort the n elements of the  array a of generic templated type Element.
/// In output the array index of type Index contains the indices of the sorted
/// array. If down is false sort in increasing order (default is decreasing
/// order).
///
/// NOTE that the array index must be created with a length >= n
/// before calling this function.
/// NOTE also that the size type for n must be the same type used for the index
/// array (templated type Index)
template <typename Element, typename Index>
void Sort(Index n, const Element* a, Index* index, Bool_t down) {
  for (Index i = 0; i < n; i++) {
    index[i] = i;
  }

  if (down) {
    std::stable_sort(index, index + n, CompareDesc<const Element*>(a));
  } else {
    std::stable_sort(index, index + n, CompareAsc<const Element*>(a));
  }
}

}  // namespace RootUtility
}  // namespace ActsExamples
