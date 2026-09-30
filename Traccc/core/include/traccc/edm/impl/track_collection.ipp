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
TRACCC_HOST_DEVICE void track<BASE>::reset_quality() {
  ndf() = {};
  chi2() = {};
  pval() = {};
  nholes() = {};
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool track<BASE>::operator==(const track<T>& other) const {
  return ((fit_outcome() == other.fit_outcome()) &&
          (params() == other.params()) && (ndf() == other.ndf()) &&
          (chi2() == other.chi2()) && (pval() == other.pval()) &&
          (nholes() == other.nholes()) &&
          (constituent_links() == other.constituent_links()));
}

}  // namespace traccc::edm
