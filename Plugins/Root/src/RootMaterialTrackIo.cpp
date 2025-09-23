// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootMaterialTrackIo.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TChain.h>
#include <TTree.h>

using namespace Acts;

void ActsPlugins::RootMaterialTrackIo::connectForRead(TChain& materialChain) {
  materialChain.SetBranchAddress("event_id", &m_eventId);
  materialChain.SetBranchAddress("v_x", &m_summaryPayload.vX);
  materialChain.SetBranchAddress("v_y", &m_summaryPayload.vY);
  materialChain.SetBranchAddress("v_z", &m_summaryPayload.vZ);
  materialChain.SetBranchAddress("v_px", &m_summaryPayload.vPx);
  materialChain.SetBranchAddress("v_py", &m_summaryPayload.vPy);
  materialChain.SetBranchAddress("v_pz", &m_summaryPayload.vPz);
  materialChain.SetBranchAddress("v_phi", &m_summaryPayload.vPhi);
  materialChain.SetBranchAddress("v_eta", &m_summaryPayload.vEta);
  materialChain.SetBranchAddress("t_X0", &m_summaryPayload.tX0);
  materialChain.SetBranchAddress("t_L0", &m_summaryPayload.tL0);
  materialChain.SetBranchAddress("mat_x", &m_stepPayload.stepXPtr);
  materialChain.SetBranchAddress("mat_y", &m_stepPayload.stepYPtr);
  materialChain.SetBranchAddress("mat_z", &m_stepPayload.stepZPtr);
  materialChain.SetBranchAddress("mat_dx", &m_stepPayload.stepDxPtr);
  materialChain.SetBranchAddress("mat_dy", &m_stepPayload.stepDyPtr);
  materialChain.SetBranchAddress("mat_dz", &m_stepPayload.stepDzPtr);
  materialChain.SetBranchAddress("mat_step_length",
                                 &m_stepPayload.stepLengthPtr);
  materialChain.SetBranchAddress("mat_X0", &m_stepPayload.stepMatX0Ptr);
  materialChain.SetBranchAddress("mat_L0", &m_stepPayload.stepMatL0Ptr);
  materialChain.SetBranchAddress("mat_A", &m_stepPayload.stepMatAPtr);
  materialChain.SetBranchAddress("mat_Z", &m_stepPayload.stepMatZPtr);
  materialChain.SetBranchAddress("mat_rho", &m_stepPayload.stepMatRhoPtr);
  if (m_cfg.surfaceInfo) {
    materialChain.SetBranchAddress("sur_id", &m_surfacePayload.surfaceIdPtr);
    materialChain.SetBranchAddress("sur_x", &m_surfacePayload.surfaceXPtr);
    materialChain.SetBranchAddress("sur_y", &m_surfacePayload.surfaceYPtr);
    materialChain.SetBranchAddress("sur_z", &m_surfacePayload.surfaceZPtr);
    materialChain.SetBranchAddress("sur_pathCorrection",
                                   &m_surfacePayload.surfacePathCorrectionPtr);
  }
}

void ActsPlugins::RootMaterialTrackIo::connectForWrite(TTree& materialTree) {
  // This sets the branch addresses for the material track
  // Set the branches
  materialTree.Branch("event_id", &m_eventId);
  materialTree.Branch("v_x", &m_summaryPayload.vX);
  materialTree.Branch("v_y", &m_summaryPayload.vY);
  materialTree.Branch("v_z", &m_summaryPayload.vZ);
  materialTree.Branch("v_px", &m_summaryPayload.vPx);
  materialTree.Branch("v_py", &m_summaryPayload.vPy);
  materialTree.Branch("v_pz", &m_summaryPayload.vPz);
  materialTree.Branch("v_phi", &m_summaryPayload.vPhi);
  materialTree.Branch("v_eta", &m_summaryPayload.vEta);
  materialTree.Branch("t_X0", &m_summaryPayload.tX0);
  materialTree.Branch("t_L0", &m_summaryPayload.tL0);
  materialTree.Branch("mat_x", &m_stepPayload.stepX);
  materialTree.Branch("mat_y", &m_stepPayload.stepY);
  materialTree.Branch("mat_z", &m_stepPayload.stepZ);
  materialTree.Branch("mat_r", &m_stepPayload.stepR);
  materialTree.Branch("mat_dx", &m_stepPayload.stepDx);
  materialTree.Branch("mat_dy", &m_stepPayload.stepDy);
  materialTree.Branch("mat_dz", &m_stepPayload.stepDz);
  materialTree.Branch("mat_step_length", &m_stepPayload.stepLength);
  materialTree.Branch("mat_X0", &m_stepPayload.stepMatX0);
  materialTree.Branch("mat_L0", &m_stepPayload.stepMatL0);
  materialTree.Branch("mat_A", &m_stepPayload.stepMatA);
  materialTree.Branch("mat_Z", &m_stepPayload.stepMatZ);
  materialTree.Branch("mat_rho", &m_stepPayload.stepMatRho);

  if (m_cfg.prePostStepInfo) {
    materialTree.Branch("mat_sx", &m_stepPayload.stepXs);
    materialTree.Branch("mat_sy", &m_stepPayload.stepYs);
    materialTree.Branch("mat_sz", &m_stepPayload.stepZs);
    materialTree.Branch("mat_ex", &m_stepPayload.stepXe);
    materialTree.Branch("mat_ey", &m_stepPayload.stepYe);
    materialTree.Branch("mat_ez", &m_stepPayload.stepZe);
  }
  if (m_cfg.surfaceInfo) {
    materialTree.Branch("sur_id", &m_surfacePayload.surfaceId);
    materialTree.Branch("sur_x", &m_surfacePayload.surfaceX);
    materialTree.Branch("sur_y", &m_surfacePayload.surfaceY);
    materialTree.Branch("sur_z", &m_surfacePayload.surfaceZ);
    materialTree.Branch("sur_r", &m_surfacePayload.surfaceR);
    materialTree.Branch("sur_distance", &m_surfacePayload.surfaceDistance);
    materialTree.Branch("sur_pathCorrection",
                        &m_surfacePayload.surfacePathCorrection);
  }
  if (m_cfg.volumeInfo) {
    materialTree.Branch("vol_id", &m_volumePayload.volumeId);
  }
}

void ActsPlugins::RootMaterialTrackIo::write(
    const GeometryContext& gctx, std::uint32_t eventNum,
    const RecordedMaterialTrack& materialTrack) {
  m_eventId = eventNum;
  // Clearing the vector first
  m_stepPayload.stepXs.clear();
  m_stepPayload.stepYs.clear();
  m_stepPayload.stepZs.clear();
  m_stepPayload.stepX.clear();
  m_stepPayload.stepY.clear();
  m_stepPayload.stepZ.clear();
  m_stepPayload.stepR.clear();
  m_stepPayload.stepXe.clear();
  m_stepPayload.stepYe.clear();
  m_stepPayload.stepZe.clear();
  m_stepPayload.stepDx.clear();
  m_stepPayload.stepDy.clear();
  m_stepPayload.stepDz.clear();
  m_stepPayload.stepLength.clear();
  m_stepPayload.stepMatX0.clear();
  m_stepPayload.stepMatL0.clear();
  m_stepPayload.stepMatA.clear();
  m_stepPayload.stepMatZ.clear();
  m_stepPayload.stepMatRho.clear();

  m_surfacePayload.surfaceId.clear();
  m_surfacePayload.surfaceX.clear();
  m_surfacePayload.surfaceY.clear();
  m_surfacePayload.surfaceZ.clear();
  m_surfacePayload.surfaceR.clear();
  m_surfacePayload.surfaceDistance.clear();
  m_surfacePayload.surfacePathCorrection.clear();

  m_volumePayload.volumeId.clear();

  auto materialInteractions = materialTrack.second.materialInteractions;

  // Reserve the vector then
  std::size_t mints = materialInteractions.size();
  m_stepPayload.stepXs.reserve(mints);
  m_stepPayload.stepYs.reserve(mints);
  m_stepPayload.stepZs.reserve(mints);
  m_stepPayload.stepX.reserve(mints);
  m_stepPayload.stepY.reserve(mints);
  m_stepPayload.stepZ.reserve(mints);
  m_stepPayload.stepR.reserve(mints);
  m_stepPayload.stepXe.reserve(mints);
  m_stepPayload.stepYe.reserve(mints);
  m_stepPayload.stepZe.reserve(mints);
  m_stepPayload.stepDx.reserve(mints);
  m_stepPayload.stepDy.reserve(mints);
  m_stepPayload.stepDz.reserve(mints);
  m_stepPayload.stepLength.reserve(mints);
  m_stepPayload.stepMatX0.reserve(mints);
  m_stepPayload.stepMatL0.reserve(mints);
  m_stepPayload.stepMatA.reserve(mints);
  m_stepPayload.stepMatZ.reserve(mints);
  m_stepPayload.stepMatRho.reserve(mints);

  m_surfacePayload.surfaceId.reserve(mints);
  m_surfacePayload.surfaceX.reserve(mints);
  m_surfacePayload.surfaceY.reserve(mints);
  m_surfacePayload.surfaceZ.reserve(mints);
  m_surfacePayload.surfaceR.reserve(mints);
  m_surfacePayload.surfaceDistance.reserve(mints);
  m_surfacePayload.surfacePathCorrection.reserve(mints);

  m_volumePayload.volumeId.reserve(mints);

  // reset the global counter
  if (m_cfg.recalculateTotals) {
    m_summaryPayload.tX0 = 0.;
    m_summaryPayload.tL0 = 0.;
  } else {
    m_summaryPayload.tX0 = materialTrack.second.materialInX0;
    m_summaryPayload.tL0 = materialTrack.second.materialInL0;
  }

  // set the track information at vertex
  m_summaryPayload.vX = materialTrack.first.first.x();
  m_summaryPayload.vY = materialTrack.first.first.y();
  m_summaryPayload.vZ = materialTrack.first.first.z();
  m_summaryPayload.vPx = materialTrack.first.second.x();
  m_summaryPayload.vPy = materialTrack.first.second.y();
  m_summaryPayload.vPz = materialTrack.first.second.z();
  m_summaryPayload.vPhi = VectorHelpers::phi(materialTrack.first.second);
  m_summaryPayload.vEta = VectorHelpers::eta(materialTrack.first.second);

  // and now loop over the material
  for (const auto& mint : materialInteractions) {
    auto direction = mint.direction.normalized();

    // The material step position information
    m_stepPayload.stepX.push_back(mint.position.x());
    m_stepPayload.stepY.push_back(mint.position.y());
    m_stepPayload.stepZ.push_back(mint.position.z());
    m_stepPayload.stepR.push_back(VectorHelpers::perp(mint.position));
    m_stepPayload.stepDx.push_back(direction.x());
    m_stepPayload.stepDy.push_back(direction.y());
    m_stepPayload.stepDz.push_back(direction.z());

    if (m_cfg.prePostStepInfo) {
      Vector3 prePos = mint.position - 0.5 * mint.pathCorrection * direction;
      Vector3 posPos = mint.position + 0.5 * mint.pathCorrection * direction;

      m_stepPayload.stepXs.push_back(prePos.x());
      m_stepPayload.stepYs.push_back(prePos.y());
      m_stepPayload.stepZs.push_back(prePos.z());
      m_stepPayload.stepXe.push_back(posPos.x());
      m_stepPayload.stepYe.push_back(posPos.y());
      m_stepPayload.stepZe.push_back(posPos.z());
    }

    // Store surface information
    if (m_cfg.surfaceInfo) {
      const Surface* surface = mint.surface;
      if (mint.intersectionID.value() != 0) {
        m_surfacePayload.surfaceId.push_back(mint.intersectionID.value());
        m_surfacePayload.surfacePathCorrection.push_back(mint.pathCorrection);
        m_surfacePayload.surfaceX.push_back(mint.intersection.x());
        m_surfacePayload.surfaceY.push_back(mint.intersection.y());
        m_surfacePayload.surfaceZ.push_back(mint.intersection.z());
        m_surfacePayload.surfaceR.push_back(
            VectorHelpers::perp(mint.intersection));
        m_surfacePayload.surfaceDistance.push_back(
            (mint.position - mint.intersection).norm());
      } else if (surface != nullptr) {
        Intersection3D sfIntersection =
            surface
                ->intersect(gctx, mint.position, mint.direction,
                            BoundaryTolerance::None())
                .closest();
        m_surfacePayload.surfaceId.push_back(surface->geometryId().value());
        m_surfacePayload.surfacePathCorrection.push_back(1.0);
        m_surfacePayload.surfaceX.push_back(sfIntersection.position().x());
        m_surfacePayload.surfaceY.push_back(sfIntersection.position().y());
        m_surfacePayload.surfaceZ.push_back(sfIntersection.position().z());
      } else {
        m_surfacePayload.surfaceId.push_back(GeometryIdentifier().value());
        m_surfacePayload.surfaceX.push_back(0);
        m_surfacePayload.surfaceY.push_back(0);
        m_surfacePayload.surfaceZ.push_back(0);
        m_surfacePayload.surfacePathCorrection.push_back(1.0);
      }
    }

    // store volume information
    if (m_cfg.volumeInfo) {
      GeometryIdentifier vlayerID;
      if (!mint.volume.empty()) {
        vlayerID = mint.volume.geometryId();
        m_volumePayload.volumeId.push_back(vlayerID.value());
      } else {
        vlayerID = vlayerID.withVolume(0)
                       .withBoundary(0)
                       .withLayer(0)
                       .withApproach(0)
                       .withSensitive(0);
        m_volumePayload.volumeId.push_back(vlayerID.value());
      }
    }

    // the material information
    const auto& mprops = mint.materialSlab;
    m_stepPayload.stepLength.push_back(mprops.thickness());
    m_stepPayload.stepMatX0.push_back(mprops.material().X0());
    m_stepPayload.stepMatL0.push_back(mprops.material().L0());
    m_stepPayload.stepMatA.push_back(mprops.material().Ar());
    m_stepPayload.stepMatZ.push_back(mprops.material().Z());
    m_stepPayload.stepMatRho.push_back(mprops.material().massDensity());
    // re-calculate if defined to do so
    if (m_cfg.recalculateTotals) {
      m_summaryPayload.tX0 += mprops.thicknessInX0();
      m_summaryPayload.tL0 += mprops.thicknessInL0();
    }
  }
}

RecordedMaterialTrack ActsPlugins::RootMaterialTrackIo::read() const {
  RecordedMaterialTrack rmTrack;
  // Fill the position and momentum
  rmTrack.first.first =
      Vector3(m_summaryPayload.vX, m_summaryPayload.vY, m_summaryPayload.vZ);
  rmTrack.first.second =
      Vector3(m_summaryPayload.vPx, m_summaryPayload.vPy, m_summaryPayload.vPz);

  // Fill the individual steps
  std::size_t msteps = m_stepPayload.stepLength.size();
  rmTrack.second.materialInteractions.reserve(msteps);
  rmTrack.second.materialInX0 = 0.;
  rmTrack.second.materialInL0 = 0.;

  for (std::size_t is = 0; is < msteps; ++is) {
    float s = m_stepPayload.stepLength[is];
    if (s == 0.) {
      continue;
    }

    float mX0 = m_stepPayload.stepMatX0[is];
    float mL0 = m_stepPayload.stepMatL0[is];

    rmTrack.second.materialInX0 += s / mX0;
    rmTrack.second.materialInL0 += s / mL0;
    /// Fill the position & the material
    MaterialInteraction mInteraction;
    mInteraction.position =
        Vector3(m_stepPayload.stepX[is], m_stepPayload.stepY[is],
                m_stepPayload.stepZ[is]);
    mInteraction.direction =
        Vector3(m_stepPayload.stepDx[is], m_stepPayload.stepDy[is],
                m_stepPayload.stepDz[is]);
    mInteraction.materialSlab = MaterialSlab(
        Material::fromMassDensity(mX0, mL0, m_stepPayload.stepMatA[is],
                                  m_stepPayload.stepMatZ[is],
                                  m_stepPayload.stepMatRho[is]),
        s);
    if (m_cfg.surfaceInfo) {
      // add the surface information to the interaction this allows the
      // mapping to be speed up
      mInteraction.intersectionID =
          GeometryIdentifier(m_surfacePayload.surfaceId[is]);
      mInteraction.intersection =
          Vector3(m_surfacePayload.surfaceX[is], m_surfacePayload.surfaceY[is],
                  m_surfacePayload.surfaceZ[is]);
      mInteraction.pathCorrection = m_surfacePayload.surfacePathCorrection[is];
    } else {
      mInteraction.intersectionID = GeometryIdentifier();
      mInteraction.intersection = Vector3(0, 0, 0);
    }
    rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
  }
  return rmTrack;
}
