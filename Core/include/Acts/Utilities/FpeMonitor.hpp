// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

class FpeMonitor {
 public:
  FpeMonitor();
  explicit FpeMonitor(int excepts);
  ~FpeMonitor();

  static void enable(int excepts);
  static void disable(int excepts);

 private:
  int m_excepts;
};

}  // namespace Acts
