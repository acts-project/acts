// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {
/**
 * Struct that sorts trackstates using their path lengths.
 * This can be used as a sorter in STL functions.
 */
struct TrackStatePathLengthSorter {
 public:
  /**
   * The sorting operator
   * @tparam identifier_t Identifier of the track state
   * @tparam parameters_t The concrete parameters type
   * @param lhs First track state
   * @param rhs Second trackstate
   * @return bool
   */
  template <typename identifier_t, typename parameters_t>
  bool operator()(const TrackState<identifier_t, parameters_t>& lhs,
                  const TrackState<identifier_t, parameters_t>& rhs) {
    return lhs.parameter.pathLength < rhs.parameter.pathLength;
  }
};
}  // namespace Acts
