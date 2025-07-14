// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootSpacePointIo.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

#include <TChain.h>
#include <TTree.h>

namespace Acts {

void RootSpacePointIo::connectForRead(
    TChain& tchain, const Experimental::SpacePointContainer2& spacePoints) {
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::X)) {
    tchain.SetBranchAddress("x", &m_x);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Y)) {
    tchain.SetBranchAddress("y", &m_y);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Z)) {
    tchain.SetBranchAddress("z", &m_z);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Time)) {
    tchain.SetBranchAddress("t", &m_t);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::R)) {
    tchain.SetBranchAddress("r", &m_r);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::VarianceZ)) {
    tchain.SetBranchAddress("var_z", &m_varZ);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::VarianceR)) {
    tchain.SetBranchAddress("var_r", &m_varR);
  }
}

void RootSpacePointIo::connectForWrite(
    TTree& ttree, const Experimental::SpacePointContainer2& spacePoints) {
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::X)) {
    ttree.Branch("x", &m_x);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Y)) {
    ttree.Branch("y", &m_y);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Z)) {
    ttree.Branch("z", &m_z);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::Time)) {
    ttree.Branch("t", &m_t);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::R)) {
    ttree.Branch("r", &m_r);
  }

  if (spacePoints.hasColumns(Experimental::SpacePointColumns::VarianceZ)) {
    ttree.Branch("var_z", &m_varZ);
  }
  if (spacePoints.hasColumns(Experimental::SpacePointColumns::VarianceR)) {
    ttree.Branch("var_r", &m_varR);
  }
}

void RootSpacePointIo::write(
    const Experimental::ConstSpacePointProxy2& spacePoint) {
  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::X)) {
    m_x = spacePoint.x();
  }
  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::Y)) {
    m_y = spacePoint.y();
  }
  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::Z)) {
    m_z = spacePoint.z();
  }

  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::Time)) {
    m_t = spacePoint.time();
  }

  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::R)) {
    m_r = spacePoint.r();
  }

  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::VarianceZ)) {
    m_varZ = spacePoint.varianceZ();
  }
  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::VarianceR)) {
    m_varR = spacePoint.varianceR();
  }
}

void RootSpacePointIo::write(
    const Experimental::SpacePointContainer2& spacePoints, TTree& ttree) {
  connectForWrite(ttree, spacePoints);

  for (Experimental::ConstSpacePointProxy2 spacePoint : spacePoints) {
    write(spacePoint);
    ttree.Fill();
  }
}

void RootSpacePointIo::read(Experimental::MutableSpacePointProxy2& spacePoint,
                            SpacePointIndex2 index) {
  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::SourceLinks)) {
    spacePoint.assignSourceLinks(std::array<SourceLink, 1>{SourceLink(index)});
  }

  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::X)) {
    spacePoint.x() = m_x;
  }
  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::Y)) {
    spacePoint.y() = m_y;
  }
  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::Z)) {
    spacePoint.z() = m_z;
  }

  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::Time)) {
    spacePoint.time() = m_t;
  }

  if (spacePoint.container().hasColumns(Experimental::SpacePointColumns::R)) {
    spacePoint.r() = m_r;
  }

  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::VarianceZ)) {
    spacePoint.varianceZ() = m_varZ;
  }
  if (spacePoint.container().hasColumns(
          Experimental::SpacePointColumns::VarianceR)) {
    spacePoint.varianceR() = m_varR;
  }
}

void RootSpacePointIo::read(TChain& tchain,
                            Experimental::SpacePointContainer2& spacePoints) {
  connectForRead(tchain, spacePoints);

  std::size_t nEntries = tchain.GetEntries();
  for (std::size_t i = 0; i < nEntries; ++i) {
    tchain.GetEntry(i);

    auto spacePoint = spacePoints.createSpacePoint();
    read(spacePoint, static_cast<SpacePointIndex2>(i));
  }
}

}  // namespace Acts
