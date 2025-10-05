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

#include "ActsExamples/Io/Root/RootMuonSpacePointWriter.hpp"

#include "TFile.h"
#include "TTree.h"
namespace ActsExamples {
RootMuonSpacePointWriter::RootMuonSpacePointWriter(const Config& config,
                                                   Acts::Logging::Level level)
    : WriterT(config.inputSpacePoints, "RootMuonSpacePointWriter", level),
      m_cfg{config} {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  m_file.reset(TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str()));
  if (m_file == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_file->cd();
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
}

RootMuonSpacePointWriter::~RootMuonSpacePointWriter() = default;
ProcessCode RootMuonSpacePointWriter::finalize() {
  m_file->cd();
  m_file->Write();
  m_file.reset();
  ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                        << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}
ProcessCode RootMuonSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const MuonSpacePointContainer& hits) {
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
