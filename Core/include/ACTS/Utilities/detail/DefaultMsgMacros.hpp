// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DEFAULT_MSG_MACROS_H
#define ACTS_DEFAULT_MSG_MACROS_H 1

// STL include(s)
#include <iostream>

namespace Acts
{
#define MSG_VERBOSE(x) // std::cout << "VERBOSE " << x << std::endl;
#define MSG_DEBUG(x)   // std::cout << "DEBUG   " << x << std::endl;
#define MSG_INFO(x)    std::cout << "INFO    " << x << std::endl;
#define MSG_ERROR(x)   std::cout << "ERROR   " << x << std::endl;
#define MSG_WARNING(x) std::cout << "WARNING " << x << std::endl;
#define MSG_FATAL(x)   std::cout << "FATAL   " << x << std::endl;
}  // end of namespace Acts

#endif // ACTS_DEFAULT_MSG_MACROS_H
