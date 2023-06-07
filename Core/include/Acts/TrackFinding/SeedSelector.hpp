// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/Seed.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

  template <typename external_spacepoint_t>
  bool voidSeedPreselector(const Acts::Seed<external_spacepoint_t>&) {
    return true;
  }

  template<typename external_spacepoint_t>
  struct SeedSelectorConfig {
    using Preselector = Delegate<bool(const Acts::Seed<external_spacepoint_t>&)>;
    Preselector preselector;

    SeedSelectorConfig();
  };

  template<typename external_spacepoint_t>
  SeedSelectorConfig<external_spacepoint_t>::SeedSelectorConfig() {
    preselector.template connect<&voidSeedPreselector<external_spacepoint_t>>();
  }
    

  template<typename external_spacepoint_t>
  class SeedSelector {
  public:
    SeedSelector(SeedSelectorConfig<external_spacepoint_t> config);
    SeedSelector(const SeedSelector&) = delete;
    SeedSelector(SeedSelector&&) = delete;
    SeedSelector& operator=(const SeedSelector&) = delete;
    SeedSelector& operator=(SeedSelector&&) = delete;
    
    ~SeedSelector() = default;

    bool passesPreSelection(const Acts::Seed<external_spacepoint_t>&) const;
    
  private:
    SeedSelectorConfig<external_spacepoint_t> m_config;
  };

  template<typename external_spacepoint_t>
  SeedSelector<external_spacepoint_t>::SeedSelector(SeedSelectorConfig<external_spacepoint_t> config)
    : m_config(std::move(config))
  {}

  template<typename external_spacepoint_t>
  bool SeedSelector<external_spacepoint_t>::passesPreSelection(const Acts::Seed<external_spacepoint_t>& seed) const
  {
    return (m_config.preselector)(seed);
  }
  
}
