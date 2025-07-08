// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootSpacePoints2Accessor.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"

#include <TChain.h>
#include <TTree.h>

namespace Acts {

void RootSpacePointAccessor::connectForRead(TChain& tchain) {
  tchain.SetBranchAddress("index", &m_index);

  tchain.SetBranchAddress("x", &m_x);
  tchain.SetBranchAddress("y", &m_y);
  tchain.SetBranchAddress("z", &m_z);

  tchain.SetBranchAddress("t", &m_t);

  tchain.SetBranchAddress("r", &m_r);

  tchain.SetBranchAddress("var_z", &m_varZ);
  tchain.SetBranchAddress("var_r", &m_varR);
}

void RootSpacePointAccessor::connectForWrite(TTree& ttree) {
  ttree.Branch("index", &m_index);

  ttree.Branch("x", &m_x);
  ttree.Branch("y", &m_y);
  ttree.Branch("z", &m_z);

  ttree.Branch("t", &m_t);

  ttree.Branch("r", &m_r);

  ttree.Branch("var_z", &m_varZ);
  ttree.Branch("var_r", &m_varR);
}

void RootSpacePointAccessor::write(
    const Experimental::ConstSpacePointProxy2& spacePoint) {
  m_index = spacePoint.index();

  m_x = spacePoint.x();
  m_y = spacePoint.y();
  m_z = spacePoint.z();

  m_t = spacePoint.time();

  m_r = spacePoint.r();

  m_varZ = spacePoint.varianceZ();
  m_varR = spacePoint.varianceR();
}

void RootSpacePointAccessor::write(
    const Experimental::SpacePointContainer2& spacePoints, TTree& ttree) {
  connectForWrite(ttree);

  for (const auto& spacePoint : spacePoints) {
    write(spacePoint);
    ttree.Fill();
  }
}

void RootSpacePointAccessor::read(
    Experimental::MutableSpacePointProxy2& spacePoint) {
  spacePoint.x() = m_x;
  spacePoint.y() = m_y;
  spacePoint.z() = m_z;

  spacePoint.time() = m_t;

  spacePoint.r() = m_r;

  spacePoint.varianceZ() = m_varZ;
  spacePoint.varianceR() = m_varR;
}

void RootSpacePointAccessor::read(
    TChain& tchain, Experimental::SpacePointContainer2& spacePoints) {
  connectForRead(tchain);

  std::size_t nEntries = tchain.GetEntries();
  for (std::size_t i = 0; i < nEntries; ++i) {
    tchain.GetEntry(i);

    auto spacePoint = spacePoints.createSpacePoint(
        std::array<SourceLink, 1>{SourceLink(i)}, 0, 0, 0);
    read(spacePoint);
  }
}

}  // namespace Acts
