// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementBase.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

/// The API has to be present though
#ifdef ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#include ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#else
#include "detail/DefaultDetectorElementBase.hpp"
#endif
