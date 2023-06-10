// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

class G4RunManager;

namespace ActsExamples {

struct Geant4Instance {
  std::mutex mutex;
  std::unique_ptr<G4RunManager> runManager;

  Geant4Instance(std::unique_ptr<G4RunManager> runManager);
  Geant4Instance(const Geant4Instance &) = delete;
  Geant4Instance &operator=(const Geant4Instance &) = delete;
  ~Geant4Instance();
};

class Geant4Manager {
 public:
  static Geant4Manager &instance();

  std::shared_ptr<Geant4Instance> create();

 private:
  Geant4Manager();
  ~Geant4Manager();

  std::weak_ptr<Geant4Instance> m_instance;
};

}  // namespace ActsExamples
