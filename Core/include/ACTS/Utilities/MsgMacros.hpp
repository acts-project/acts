// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MSG_MACROS_H
#define ACTS_MSG_MACROS_H 1

#ifdef ACTS_MSG_MACROS_PLUGIN
#include ACTS_MSG_MACROS_PLUGIN
#else
static_assert(false, "no message macros defined");
#endif

#endif  // ACTS_MSG_MACROS_H
