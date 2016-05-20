// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinStepping.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINSTEPPING_H
#define ACTS_GEOMETRYUTILS_BINSTEPPING_H 1

/** Very simple binStepping class to deal with open and closed bins */

namespace Acts {
    
    static bool decrement(size_t& bin, size_t nBins, bool openBin) {
        if (openBin && !bin) return false;
        bin = bin ? bin-1 : nBins-1;
        return true; 
    }
    
    static bool increment(size_t& bin, size_t nBins, bool openBin) {
        if (openBin && (bin+1==nBins)) return false;
        bin = (bin==nBins-1) ? 0 : bin+1;
        return true; 
    }
}

#endif
