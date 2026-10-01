// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Mark this as a "system header". To suppress all warnings from oneDPL.
// This is needed because at the time of writing we cannot provide oneDPL with
// "-isystem" to the oneAPI compiler.
#pragma clang system_header

// oneDPL include(s).
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/numeric>
