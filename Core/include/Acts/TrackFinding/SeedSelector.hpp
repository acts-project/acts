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
  bool voidSeedSelector(const Acts::Seed<external_spacepoint_t>&) {
    // By default always return true -> no selection
    return true;
  }

  // Seed Selector options
  template<typename external_spacepoint_t>
  struct SeedSelectorOptions {
    using Selector = Delegate<bool(const Acts::Seed<external_spacepoint_t>&)>;
    Selector selector;

    SeedSelectorOptions();
  };

  template<typename external_spacepoint_t>
  SeedSelectorOptions<external_spacepoint_t>::SeedSelectorOptions() {
    selector.template connect<&voidSeedSelector<external_spacepoint_t>>();
  }
    

  // Seed Selector
  template<typename external_spacepoint_t>
  class SeedSelector {
  public:
    SeedSelector(SeedSelectorOptions<external_spacepoint_t>&& options);
    SeedSelector(const SeedSelector&) = delete;
    SeedSelector(SeedSelector&&) = delete;
    SeedSelector& operator=(const SeedSelector&) = delete;
    SeedSelector& operator=(SeedSelector&&) = delete;
    
    ~SeedSelector() = default;

    bool passesQualitySelection(const Acts::Seed<external_spacepoint_t>&) const;
    
  private:
    SeedSelectorOptions<external_spacepoint_t> m_options;
  };

  template<typename external_spacepoint_t>
  SeedSelector<external_spacepoint_t>::SeedSelector(SeedSelectorOptions<external_spacepoint_t>&& options)
    : m_options(std::move(options))
  {}

  template<typename external_spacepoint_t>
  bool SeedSelector<external_spacepoint_t>::passesQualitySelection(const Acts::Seed<external_spacepoint_t>& seed) const
  {
    return (m_options.selector)(seed);
  }
  
}
