// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc {

/// Enumeration of strategies to await for completion of asynchronous
/// operations.
enum class await_strategy {
  sync_event,  ///< Synchronous waiting for an event to complete
  callback     ///< Suspension with a callback function
};

}  // namespace traccc
