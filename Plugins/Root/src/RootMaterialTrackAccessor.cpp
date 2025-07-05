// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMaterialTrackAccessor.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TChain.h>
#include <TTree.h>

void Acts::RootMaterialTrackAccessor::connectForRead(TChain& materialChain) {
  materialChain.SetBranchAddress("event_id", &m_eventId);
  materialChain.SetBranchAddress("v_x", &m_vX);
  materialChain.SetBranchAddress("v_y", &m_vY);
  materialChain.SetBranchAddress("v_z", &m_vZ);
  materialChain.SetBranchAddress("v_px", &m_vPx);
  materialChain.SetBranchAddress("v_py", &m_vPy);
  materialChain.SetBranchAddress("v_pz", &m_vPz);
  materialChain.SetBranchAddress("v_phi", &m_vPhi);
  materialChain.SetBranchAddress("v_eta", &m_vEta);
  materialChain.SetBranchAddress("t_X0", &m_tX0);
  materialChain.SetBranchAddress("t_L0", &m_tL0);
  materialChain.SetBranchAddress("mat_x", &m_stepXPtr);
  materialChain.SetBranchAddress("mat_y", &m_stepYPtr);
  materialChain.SetBranchAddress("mat_z", &m_stepZPtr);
  materialChain.SetBranchAddress("mat_dx", &m_stepDxPtr);
  materialChain.SetBranchAddress("mat_dy", &m_stepDyPtr);
  materialChain.SetBranchAddress("mat_dz", &m_stepDzPtr);
  materialChain.SetBranchAddress("mat_step_length", &m_stepLengthPtr);
  materialChain.SetBranchAddress("mat_X0", &m_stepMatX0Ptr);
  materialChain.SetBranchAddress("mat_L0", &m_stepMatL0Ptr);
  materialChain.SetBranchAddress("mat_A", &m_stepMatAPtr);
  materialChain.SetBranchAddress("mat_Z", &m_stepMatZPtr);
  materialChain.SetBranchAddress("mat_rho", &m_stepMatRhoPtr);
  if (m_cfg.surfaceInfo) {
    materialChain.SetBranchAddress("sur_id", &m_surfaceIdPtr);
    materialChain.SetBranchAddress("sur_x", &m_surfaceXPtr);
    materialChain.SetBranchAddress("sur_y", &m_surfaceYPtr);
    materialChain.SetBranchAddress("sur_z", &m_surfaceZPtr);
    materialChain.SetBranchAddress("sur_pathCorrection",
                                   &m_surfacePathCorrectionPtr);
  }
}

void Acts::RootMaterialTrackAccessor::connectForWrite(TTree& materialTree) {
  // This sets the branch addresses for the material track
  // Set the branches
  materialTree.Branch("event_id", &m_eventId);
  materialTree.Branch("v_x", &m_vX);
  materialTree.Branch("v_y", &m_vY);
  materialTree.Branch("v_z", &m_vZ);
  materialTree.Branch("v_px", &m_vPx);
  materialTree.Branch("v_py", &m_vPy);
  materialTree.Branch("v_pz", &m_vPz);
  materialTree.Branch("v_phi", &m_vPhi);
  materialTree.Branch("v_eta", &m_vEta);
  materialTree.Branch("t_X0", &m_tX0);
  materialTree.Branch("t_L0", &m_tL0);
  materialTree.Branch("mat_x", &m_stepX);
  materialTree.Branch("mat_y", &m_stepY);
  materialTree.Branch("mat_z", &m_stepZ);
  materialTree.Branch("mat_r", &m_stepR);
  materialTree.Branch("mat_dx", &m_stepDx);
  materialTree.Branch("mat_dy", &m_stepDy);
  materialTree.Branch("mat_dz", &m_stepDz);
  materialTree.Branch("mat_step_length", &m_stepLength);
  materialTree.Branch("mat_X0", &m_stepMatX0);
  materialTree.Branch("mat_L0", &m_stepMatL0);
  materialTree.Branch("mat_A", &m_stepMatA);
  materialTree.Branch("mat_Z", &m_stepMatZ);
  materialTree.Branch("mat_rho", &m_stepMatRho);

  if (m_cfg.prePostStepInfo) {
    materialTree.Branch("mat_sx", &m_stepXs);
    materialTree.Branch("mat_sy", &m_stepYs);
    materialTree.Branch("mat_sz", &m_stepZs);
    materialTree.Branch("mat_ex", &m_stepXe);
    materialTree.Branch("mat_ey", &m_stepYe);
    materialTree.Branch("mat_ez", &m_stepZe);
  }
  if (m_cfg.surfaceInfo) {
    materialTree.Branch("sur_id", &m_surfaceId);
    materialTree.Branch("sur_x", &m_surfaceX);
    materialTree.Branch("sur_y", &m_surfaceY);
    materialTree.Branch("sur_z", &m_surfaceZ);
    materialTree.Branch("sur_r", &m_surfaceR);
    materialTree.Branch("sur_distance", &m_surfaceDistance);
    materialTree.Branch("sur_pathCorrection", &m_surfacePathCorrection);
  }
  if (m_cfg.volumeInfo) {
    materialTree.Branch("vol_id", &m_volumeId);
  }
}

void Acts::RootMaterialTrackAccessor::write(
    const GeometryContext& gctx, std::uint32_t eventNum,
    const RecordedMaterialTrack& materialTrack) {
  m_eventId = eventNum;
  // Clearing the vector first
  m_stepXs.clear();
  m_stepYs.clear();
  m_stepZs.clear();
  m_stepX.clear();
  m_stepY.clear();
  m_stepZ.clear();
  m_stepR.clear();
  m_stepXe.clear();
  m_stepYe.clear();
  m_stepZe.clear();
  m_stepDx.clear();
  m_stepDy.clear();
  m_stepDz.clear();
  m_stepLength.clear();
  m_stepMatX0.clear();
  m_stepMatL0.clear();
  m_stepMatA.clear();
  m_stepMatZ.clear();
  m_stepMatRho.clear();

  m_surfaceId.clear();
  m_surfaceX.clear();
  m_surfaceY.clear();
  m_surfaceZ.clear();
  m_surfaceR.clear();
  m_surfaceDistance.clear();
  m_surfacePathCorrection.clear();

  m_volumeId.clear();

  auto materialInteractions = materialTrack.second.materialInteractions;

  // Reserve the vector then
  std::size_t mints = materialInteractions.size();
  m_stepXs.reserve(mints);
  m_stepYs.reserve(mints);
  m_stepZs.reserve(mints);
  m_stepX.reserve(mints);
  m_stepY.reserve(mints);
  m_stepZ.reserve(mints);
  m_stepR.reserve(mints);
  m_stepXe.reserve(mints);
  m_stepYe.reserve(mints);
  m_stepZe.reserve(mints);
  m_stepDx.reserve(mints);
  m_stepDy.reserve(mints);
  m_stepDz.reserve(mints);
  m_stepLength.reserve(mints);
  m_stepMatX0.reserve(mints);
  m_stepMatL0.reserve(mints);
  m_stepMatA.reserve(mints);
  m_stepMatZ.reserve(mints);
  m_stepMatRho.reserve(mints);

  m_surfaceId.reserve(mints);
  m_surfaceX.reserve(mints);
  m_surfaceY.reserve(mints);
  m_surfaceZ.reserve(mints);
  m_surfaceR.reserve(mints);
  m_surfaceDistance.reserve(mints);
  m_surfacePathCorrection.reserve(mints);

  m_volumeId.reserve(mints);

  // reset the global counter
  if (m_cfg.recalculateTotals) {
    m_tX0 = 0.;
    m_tL0 = 0.;
  } else {
    m_tX0 = materialTrack.second.materialInX0;
    m_tL0 = materialTrack.second.materialInL0;
  }

  // set the track information at vertex
  m_vX = materialTrack.first.first.x();
  m_vY = materialTrack.first.first.y();
  m_vZ = materialTrack.first.first.z();
  m_vPx = materialTrack.first.second.x();
  m_vPy = materialTrack.first.second.y();
  m_vPz = materialTrack.first.second.z();
  m_vPhi = VectorHelpers::phi(materialTrack.first.second);
  m_vEta = VectorHelpers::eta(materialTrack.first.second);

  // and now loop over the material
  for (const auto& mint : materialInteractions) {
    auto direction = mint.direction.normalized();

    // The material step position information
    m_stepX.push_back(mint.position.x());
    m_stepY.push_back(mint.position.y());
    m_stepZ.push_back(mint.position.z());
    m_stepR.push_back(VectorHelpers::perp(mint.position));
    m_stepDx.push_back(direction.x());
    m_stepDy.push_back(direction.y());
    m_stepDz.push_back(direction.z());

    if (m_cfg.prePostStepInfo) {
      Acts::Vector3 prePos =
          mint.position - 0.5 * mint.pathCorrection * direction;
      Acts::Vector3 posPos =
          mint.position + 0.5 * mint.pathCorrection * direction;

      m_stepXs.push_back(prePos.x());
      m_stepYs.push_back(prePos.y());
      m_stepZs.push_back(prePos.z());
      m_stepXe.push_back(posPos.x());
      m_stepYe.push_back(posPos.y());
      m_stepZe.push_back(posPos.z());
    }

    // Store surface information
    if (m_cfg.surfaceInfo) {
      const Acts::Surface* surface = mint.surface;
      if (mint.intersectionID.value() != 0) {
        m_surfaceId.push_back(mint.intersectionID.value());
        m_surfacePathCorrection.push_back(mint.pathCorrection);
        m_surfaceX.push_back(mint.intersection.x());
        m_surfaceY.push_back(mint.intersection.y());
        m_surfaceZ.push_back(mint.intersection.z());
        m_surfaceR.push_back(VectorHelpers::perp(mint.intersection));
        m_surfaceDistance.push_back((mint.position - mint.intersection).norm());
      } else if (surface != nullptr) {
        auto sfIntersection =
            surface
                ->intersect(gctx, mint.position, mint.direction,
                            Acts::BoundaryTolerance::None())
                .closest();
        m_surfaceId.push_back(surface->geometryId().value());
        m_surfacePathCorrection.push_back(1.0);
        m_surfaceX.push_back(sfIntersection.position().x());
        m_surfaceY.push_back(sfIntersection.position().y());
        m_surfaceZ.push_back(sfIntersection.position().z());
      } else {
        m_surfaceId.push_back(Acts::GeometryIdentifier().value());
        m_surfaceX.push_back(0);
        m_surfaceY.push_back(0);
        m_surfaceZ.push_back(0);
        m_surfacePathCorrection.push_back(1.0);
      }
    }

    // store volume information
    if (m_cfg.volumeInfo) {
      Acts::GeometryIdentifier vlayerID;
      if (!mint.volume.empty()) {
        vlayerID = mint.volume.geometryId();
        m_volumeId.push_back(vlayerID.value());
      } else {
        vlayerID = vlayerID.withVolume(0)
                       .withBoundary(0)
                       .withLayer(0)
                       .withApproach(0)
                       .withSensitive(0);
        m_volumeId.push_back(vlayerID.value());
      }
    }

    // the material information
    const auto& mprops = mint.materialSlab;
    m_stepLength.push_back(mprops.thickness());
    m_stepMatX0.push_back(mprops.material().X0());
    m_stepMatL0.push_back(mprops.material().L0());
    m_stepMatA.push_back(mprops.material().Ar());
    m_stepMatZ.push_back(mprops.material().Z());
    m_stepMatRho.push_back(mprops.material().massDensity());
    // re-calculate if defined to do so
    if (m_cfg.recalculateTotals) {
      m_tX0 += mprops.thicknessInX0();
      m_tL0 += mprops.thicknessInL0();
    }
  }
}

Acts::RecordedMaterialTrack Acts::RootMaterialTrackAccessor::read() const {
  Acts::RecordedMaterialTrack rmTrack;
  // Fill the position and momentum
  rmTrack.first.first = Acts::Vector3(m_vX, m_vY, m_vZ);
  rmTrack.first.second = Acts::Vector3(m_vPx, m_vPy, m_vPz);

  // Fill the individual steps
  std::size_t msteps = m_stepLength.size();
  rmTrack.second.materialInteractions.reserve(msteps);
  rmTrack.second.materialInX0 = 0.;
  rmTrack.second.materialInL0 = 0.;

  for (std::size_t is = 0; is < msteps; ++is) {
    float s = m_stepLength[is];
    if (s == 0.) {
      continue;
    }

    float mX0 = m_stepMatX0[is];
    float mL0 = m_stepMatL0[is];

    rmTrack.second.materialInX0 += s / mX0;
    rmTrack.second.materialInL0 += s / mL0;
    /// Fill the position & the material
    Acts::MaterialInteraction mInteraction;
    mInteraction.position =
        Acts::Vector3(m_stepX[is], m_stepY[is], m_stepZ[is]);
    mInteraction.direction =
        Acts::Vector3(m_stepDx[is], m_stepDy[is], m_stepDz[is]);
    mInteraction.materialSlab = Acts::MaterialSlab(
        Acts::Material::fromMassDensity(mX0, mL0, m_stepMatA[is],
                                        m_stepMatZ[is], m_stepMatRho[is]),
        s);
    if (m_cfg.surfaceInfo) {
      // add the surface information to the interaction this allows the
      // mapping to be speed up
      mInteraction.intersectionID = Acts::GeometryIdentifier(m_surfaceId[is]);
      mInteraction.intersection =
          Acts::Vector3(m_surfaceX[is], m_surfaceY[is], m_surfaceZ[is]);
      mInteraction.pathCorrection = m_surfacePathCorrection[is];
    } else {
      mInteraction.intersectionID = Acts::GeometryIdentifier();
      mInteraction.intersection = Acts::Vector3(0, 0, 0);
    }
    rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
  }
  return rmTrack;
}
