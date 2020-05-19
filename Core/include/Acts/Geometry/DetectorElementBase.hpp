// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/// This is the plugin mechanism to exchange the entire DetectorElementBase
///
/// By defining ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT pre-compile time the
/// detector element entire detector element can be exchanged with a file
/// provided by the client.
///
/// The API has to be present though
#ifdef ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#include ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#else
#include "detail/DefaultDetectorElementBase.hpp"
#endif
