// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Plugins/Root/RootSpacePoints2Accessor.hpp"

#include <TChain.h>
#include <TTree.h>

namespace Acts {

void RootSpacePointIo::connectForRead(
    TChain& tchain, const Experimental::SpacePointContainer2& spacePoints) {
  tchain.SetBranchAddress("index", &m_index);

  tchain.SetBranchAddress("x", &m_x);
  tchain.SetBranchAddress("y", &m_y);
  tchain.SetBranchAddress("z", &m_z);

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::Time)) {
    tchain.SetBranchAddress("t", &m_t);
  }

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::R)) {
    tchain.SetBranchAddress("r", &m_r);
  }

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceZ)) {
    tchain.SetBranchAddress("var_z", &m_varZ);
  }
  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceR)) {
    tchain.SetBranchAddress("var_r", &m_varR);
  }
}

void RootSpacePointIo::connectForWrite(
    TTree& ttree, const Experimental::SpacePointContainer2& spacePoints) {
  ttree.Branch("index", &m_index);

  ttree.Branch("x", &m_x);
  ttree.Branch("y", &m_y);
  ttree.Branch("z", &m_z);

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::Time)) {
    ttree.Branch("t", &m_t);
  }

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::R)) {
    ttree.Branch("r", &m_r);
  }

  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceZ)) {
    ttree.Branch("var_z", &m_varZ);
  }
  if (spacePoints.hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceR)) {
    ttree.Branch("var_r", &m_varR);
  }
}

void RootSpacePointIo::write(
    const Experimental::ConstSpacePointProxy2& spacePoint) {
  m_index = spacePoint.index();

  m_x = spacePoint.x();
  m_y = spacePoint.y();
  m_z = spacePoint.z();

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::Time)) {
    m_t = spacePoint.time();
  }

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::R)) {
    m_r = spacePoint.r();
  }

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceZ)) {
    m_varZ = spacePoint.varianceZ();
  }
  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceR)) {
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

void RootSpacePointIo::read(Experimental::MutableSpacePointProxy2& spacePoint) {
  spacePoint.x() = m_x;
  spacePoint.y() = m_y;
  spacePoint.z() = m_z;

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::Time)) {
    spacePoint.time() = m_t;
  }

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::R)) {
    spacePoint.r() = m_r;
  }

  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceZ)) {
    spacePoint.varianceZ() = m_varZ;
  }
  if (spacePoint.container().hasExtraColumns(
          Experimental::SpacePointKnownExtraColumn::VarianceR)) {
    spacePoint.varianceR() = m_varR;
  }
}

void RootSpacePointIo::read(TChain& tchain,
                            Experimental::SpacePointContainer2& spacePoints) {
  connectForRead(tchain, spacePoints);

  std::size_t nEntries = tchain.GetEntries();
  for (std::size_t i = 0; i < nEntries; ++i) {
    tchain.GetEntry(i);

    auto spacePoint = spacePoints.createSpacePoint(
        std::array<SourceLink, 1>{SourceLink(i)}, 0, 0, 0);
    read(spacePoint);
  }
}

}  // namespace Acts
