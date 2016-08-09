// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinningType.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNINGTYPE_H
#define ACTS_GEOMETRYUTILS_BINNINGTYPE_H 1

// STL include(s)
#include <string>
#include <vector>

namespace Acts {

/** @enum BinningType, BinningOption & BinningAccess

   - BinningType:

     Enumeration to qualify the binning type for the use of the
     LayerArrayCreator and the TrackingVolumeArrayCreator

    - BinningOption:
      open:   [0,max]
      closed:  0 -> nextbin -> max -> 0

    - BinningValue
      necessary access to global positions

   */
enum BinningType { equidistant, arbitrary };

/** enum BinValue */
enum BinningOption { open, closed };

/**  how to take the global / local position */
enum BinningValue { binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta };

/** screen output option */
static std::vector<std::string> binningValueNames
    = {"binX", "binY", "binZ", "binR", "binPhi", "binRPhi", "binH", "binEta"};
}
#endif  // ACTS_GEOMETRYUTILS_BINNINGTYPE_H
