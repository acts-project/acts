// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESDISCBOUNDS_H
#define ACTS_SURFACESDISCBOUNDS_H

#include "ACTS/Surfaces/SurfaceBounds.hpp"

namespace Acts {

/// @class DiscBounds
///
/// common base class for all bounds that are in a r/phi frame
///  - simply introduced to avoid wrong bound assigments to surfaces

class DiscBounds : public SurfaceBounds
{
public:
  /// Default Constructor
  ///
  /// @param sSize is the size of the store
  DiscBounds(size_t sSize = 0) : SurfaceBounds(sSize) {}

  /// Destructor
  virtual ~DiscBounds() {}

  /// Virtual Constructor
  virtual DiscBounds*
  clone() const = 0;
};

}  // end of namespace

#endif  // ACTS_SURFACESDISCBOUNDS_H
