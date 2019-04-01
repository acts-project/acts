// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinningType.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <string>
#include <vector>

namespace Acts {

/// @enum BinningType, BinningOption & BinningAccess
///
///- BinningType:
///
///  Enumeration to qualify the binning type for the use of the
///  LayerArrayCreator and the TrackingVolumeArrayCreator
///
/// - BinningOption:
///   open:   [0,max]
///   closed:  0 -> nextbin -> max -> 0
///
/// - BinningValue
///   necessary access to global positions
///
enum BinningType { equidistant, arbitrary };

/// @brief flag for open/closed bins
enum BinningOption { open, closed };

/// @enum BinningValue how to take the global / local position
enum BinningValue {
  binX,
  binY,
  binZ,
  binR,
  binPhi,
  binRPhi,
  binH,
  binEta,
  binMag
};

/// @brief screen output option
static const std::vector<std::string> binningValueNames = {"binX",
                                                           "binY",
                                                           "binZ",
                                                           "binR",
                                                           "binPhi",
                                                           "binRPhi",
                                                           "binH",
                                                           "binEta",
                                                           "binMag"};
}