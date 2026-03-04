// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootSpacePointIo.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

#include <TChain.h>
#include <TTree.h>

using namespace Acts;

namespace ActsPlugins {

void RootSpacePointIo::connectForRead(TChain& tchain,
                                      const SpacePointContainer2& spacePoints) {
  using enum SpacePointColumns;

  if (spacePoints.hasColumns(X)) {
    tchain.SetBranchAddress("x", &m_x);
  }
  if (spacePoints.hasColumns(Y)) {
    tchain.SetBranchAddress("y", &m_y);
  }
  if (spacePoints.hasColumns(Z)) {
    tchain.SetBranchAddress("z", &m_z);
  }

  if (spacePoints.hasColumns(Time)) {
    tchain.SetBranchAddress("t", &m_t);
  }

  if (spacePoints.hasColumns(R)) {
    tchain.SetBranchAddress("r", &m_r);
  }

  if (spacePoints.hasColumns(VarianceZ)) {
    tchain.SetBranchAddress("var_z", &m_varZ);
  }
  if (spacePoints.hasColumns(VarianceR)) {
    tchain.SetBranchAddress("var_r", &m_varR);
  }
}

void RootSpacePointIo::connectForWrite(
    TTree& ttree, const SpacePointContainer2& spacePoints) {
  using enum SpacePointColumns;

  if (spacePoints.hasColumns(X)) {
    ttree.Branch("x", &m_x);
  }
  if (spacePoints.hasColumns(Y)) {
    ttree.Branch("y", &m_y);
  }
  if (spacePoints.hasColumns(Z)) {
    ttree.Branch("z", &m_z);
  }

  if (spacePoints.hasColumns(Time)) {
    ttree.Branch("t", &m_t);
  }

  if (spacePoints.hasColumns(R)) {
    ttree.Branch("r", &m_r);
  }

  if (spacePoints.hasColumns(VarianceZ)) {
    ttree.Branch("var_z", &m_varZ);
  }
  if (spacePoints.hasColumns(VarianceR)) {
    ttree.Branch("var_r", &m_varR);
  }
}

void RootSpacePointIo::write(const ConstSpacePointProxy2& spacePoint) {
  using enum SpacePointColumns;

  if (spacePoint.container().hasColumns(X)) {
    m_x = spacePoint.x();
  }
  if (spacePoint.container().hasColumns(Y)) {
    m_y = spacePoint.y();
  }
  if (spacePoint.container().hasColumns(Z)) {
    m_z = spacePoint.z();
  }

  if (spacePoint.container().hasColumns(Time)) {
    m_t = spacePoint.time();
  }

  if (spacePoint.container().hasColumns(R)) {
    m_r = spacePoint.r();
  }

  if (spacePoint.container().hasColumns(VarianceZ)) {
    m_varZ = spacePoint.varianceZ();
  }
  if (spacePoint.container().hasColumns(VarianceR)) {
    m_varR = spacePoint.varianceR();
  }
}

void RootSpacePointIo::write(const SpacePointContainer2& spacePoints,
                             TTree& ttree) {
  connectForWrite(ttree, spacePoints);

  for (ConstSpacePointProxy2 spacePoint : spacePoints) {
    write(spacePoint);
    ttree.Fill();
  }
}

void RootSpacePointIo::read(MutableSpacePointProxy2& spacePoint,
                            SpacePointIndex2 index) {
  using enum SpacePointColumns;

  if (spacePoint.container().hasColumns(SourceLinks)) {
    spacePoint.assignSourceLinks(std::array<SourceLink, 1>{SourceLink(index)});
  }

  if (spacePoint.container().hasColumns(X)) {
    spacePoint.x() = m_x;
  }
  if (spacePoint.container().hasColumns(Y)) {
    spacePoint.y() = m_y;
  }
  if (spacePoint.container().hasColumns(Z)) {
    spacePoint.z() = m_z;
  }

  if (spacePoint.container().hasColumns(Time)) {
    spacePoint.time() = m_t;
  }

  if (spacePoint.container().hasColumns(R)) {
    spacePoint.r() = m_r;
  }

  if (spacePoint.container().hasColumns(VarianceZ)) {
    spacePoint.varianceZ() = m_varZ;
  }
  if (spacePoint.container().hasColumns(VarianceR)) {
    spacePoint.varianceR() = m_varR;
  }
}

void RootSpacePointIo::read(TChain& tchain, SpacePointContainer2& spacePoints) {
  connectForRead(tchain, spacePoints);

  std::size_t nEntries = tchain.GetEntries();
  for (std::size_t i = 0; i < nEntries; ++i) {
    tchain.GetEntry(i);

    auto spacePoint = spacePoints.createSpacePoint();
    read(spacePoint, static_cast<SpacePointIndex2>(i));
  }
}

}  // namespace ActsPlugins
