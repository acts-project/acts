// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>

class G4RunManager;
class G4VUserPhysicsList;

namespace ActsExamples {

class PhysicsListFactory;

struct Geant4Instance {
  std::mutex mutex;
  std::unique_ptr<G4RunManager> runManager;
  std::unordered_map<std::string, G4VUserPhysicsList *> physicsLists;
  G4VUserPhysicsList *anonymousPhysicList{};

  Geant4Instance(std::unique_ptr<G4RunManager> runManager);
  Geant4Instance(const Geant4Instance &) = delete;
  Geant4Instance &operator=(const Geant4Instance &) = delete;
  ~Geant4Instance();

  void registerPhysicsList(std::string name,
                           std::unique_ptr<G4VUserPhysicsList> physicsList);
  G4VUserPhysicsList *getPhysicsList(const std::string &name) const;
  G4VUserPhysicsList *createRegisterAndGetPhysicsList(const std::string &name);
  G4VUserPhysicsList *registerAndGetAnonymousPhysicsList(
      std::unique_ptr<G4VUserPhysicsList> physicsList);

  void tweekLogging(int level) const;

  void clearPhysicsList();
};

class Geant4Manager {
 public:
  static Geant4Manager &instance();

  /// This can only be called once due to Geant4 limitations
  std::shared_ptr<Geant4Instance> create();

  void registerPhysicsListFactory(
      std::string name, std::shared_ptr<PhysicsListFactory> physicsListFactroy);
  std::unique_ptr<G4VUserPhysicsList> createPhysicsList(
      const std::string &name) const;

 private:
  Geant4Manager();
  ~Geant4Manager();

  bool m_created = false;
  std::weak_ptr<Geant4Instance> m_instance;
  std::unordered_map<std::string, std::shared_ptr<PhysicsListFactory>>
      m_physicsListsFactory;
};

}  // namespace ActsExamples
