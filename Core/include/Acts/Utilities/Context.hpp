// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/any.hpp>

namespace Acts {

/// This is the central definition of the Acts
/// payload object that is propagated through the code
/// to allow for event/thread context

using Context        = boost::any;
using DefaultContext = boost::any;

}  // namespace Acts
