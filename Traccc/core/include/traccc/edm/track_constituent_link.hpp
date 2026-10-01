// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <ostream>

namespace traccc::edm {

/// A link to a constituent of a track
///
/// It is a poor man's version of @c Acts::SourceLink. Only providing the
/// functionality needed in the GPU algorithms.
///
struct track_constituent_link {
  /// The type of the constituent
  enum constituent_type : unsigned short {
    measurement = 0,  ///< The link points at a measurement
    track_state = 1   ///< The link points at a track state
  };

  unsigned short type;  ///< The type of the constituent
  unsigned int index;   ///< The index of the constituent in its collection

  /// Equality operator
  ///
  /// For some reason Clang fails to generate it automatically for this type.
  ///
  bool operator==(const track_constituent_link&) const = default;

 private:
  /// @returns a string stream that prints the constituent link
  TRACCC_HOST
  friend std::ostream& operator<<(std::ostream& os,
                                  const track_constituent_link& l) {
    if (l.type == constituent_type::measurement) {
      os << "measurement index: " << l.index;
    } else if (l.type == constituent_type::track_state) {
      os << "track state index: " << l.index;
    } else {
      os << "ERROR: Unknown link type!";
    }

    return os;
  }

};  // struct track_constituent_link

}  // namespace traccc::edm
