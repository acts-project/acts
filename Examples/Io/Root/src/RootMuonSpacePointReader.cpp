// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMuonSpacePointReader.hpp"


using namespace Acts;
using namespace Acts::UnitLiterals;
namespace ActsExamples {

   RootMuonSpacePointReader::RootMuonSpacePointReader(const Config& config, Acts::Logging::Level level):
    m_cfg{config, getDefaultLogger("RootMuonSpacePointReader", level)}{

    }

  otMuonSpacePointReader::  ~RootMuonSpacePointReader() =default;

  std::pair<std::size_t, std::size_t> RootMuonSpacePointReader::availableEvents() const {
      return std::make_pair(0u,0u);
  }

  ProcessCode RootMuonSpacePointReader::read(const ActsExamples::AlgorithmContext &context) {
     return ProcessCode::SUCCESS;
  }


  
}  // namespace ActsExamples
