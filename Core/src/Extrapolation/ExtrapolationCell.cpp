// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationCell.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Extrapolation/ExtrapolationCell.hpp"

std::vector<std::string> Acts::ExtrapolationCode::s_ecodeNames
    = {"Unset",
       "InProgress",
       "SuccessDestination",
       "SuccessBoundaryReached",
       "SuccessPathLimit",
       "SuccessMaterialLimit",
       "Recovered",
       "FailureDestination",
       "FailureLoop",
       "FailureNavigation",
       "FailureUpdateKill",
       "FailureConfiguration",
       "LeftKnownWorld"};
