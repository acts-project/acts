// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMaterialTrackPayload.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TChain.h>
#include <TTree.h>

Acts::RootMaterialTrackPayload::~RootMaterialTrackPayload() {
  delete m_step_sx;
  delete m_step_sy;
  delete m_step_sz;
  delete m_step_ex;
  delete m_step_ey;
  delete m_step_ez;
  delete m_step_x;
  delete m_step_y;
  delete m_step_z;
  delete m_step_r;
  delete m_step_dx;
  delete m_step_dy;
  delete m_step_dz;
  delete m_step_length;
  delete m_step_X0;
  delete m_step_L0;
  delete m_step_A;
  delete m_step_Z;
  delete m_step_rho;
  delete m_sur_id;
  delete m_sur_x;
  delete m_sur_y;
  delete m_sur_z;
  delete m_sur_r;
  delete m_sur_distance;
  delete m_sur_pathCorrection;
  delete m_vol_id;
}

void Acts::RootMaterialTrackPayload::connectForRead(TChain& materialChain) {
  materialChain.SetBranchAddress("event_id", &m_eventId);
  materialChain.SetBranchAddress("v_x", &m_v_x);
  materialChain.SetBranchAddress("v_y", &m_v_y);
  materialChain.SetBranchAddress("v_z", &m_v_z);
  materialChain.SetBranchAddress("v_px", &m_v_px);
  materialChain.SetBranchAddress("v_py", &m_v_py);
  materialChain.SetBranchAddress("v_pz", &m_v_pz);
  materialChain.SetBranchAddress("v_phi", &m_v_phi);
  materialChain.SetBranchAddress("v_eta", &m_v_eta);
  materialChain.SetBranchAddress("t_X0", &m_tX0);
  materialChain.SetBranchAddress("t_L0", &m_tL0);
  materialChain.SetBranchAddress("mat_x", &m_step_x);
  materialChain.SetBranchAddress("mat_y", &m_step_y);
  materialChain.SetBranchAddress("mat_z", &m_step_z);
  materialChain.SetBranchAddress("mat_dx", &m_step_dx);
  materialChain.SetBranchAddress("mat_dy", &m_step_dy);
  materialChain.SetBranchAddress("mat_dz", &m_step_dz);
  materialChain.SetBranchAddress("mat_step_length", &m_step_length);
  materialChain.SetBranchAddress("mat_X0", &m_step_X0);
  materialChain.SetBranchAddress("mat_L0", &m_step_L0);
  materialChain.SetBranchAddress("mat_A", &m_step_A);
  materialChain.SetBranchAddress("mat_Z", &m_step_Z);
  materialChain.SetBranchAddress("mat_rho", &m_step_rho);
  if (m_surfaceInfo) {
    materialChain.SetBranchAddress("sur_id", &m_sur_id);
    materialChain.SetBranchAddress("sur_x", &m_sur_x);
    materialChain.SetBranchAddress("sur_y", &m_sur_y);
    materialChain.SetBranchAddress("sur_z", &m_sur_z);
    materialChain.SetBranchAddress("sur_pathCorrection", &m_sur_pathCorrection);
  }
}

void Acts::RootMaterialTrackPayload::connectForWrite(TTree& materialTree) {
  // This sets the branch addresses for the material track
  // Set the branches
  materialTree.Branch("event_id", &m_eventId);
  materialTree.Branch("v_x", &m_v_x);
  materialTree.Branch("v_y", &m_v_y);
  materialTree.Branch("v_z", &m_v_z);
  materialTree.Branch("v_px", &m_v_px);
  materialTree.Branch("v_py", &m_v_py);
  materialTree.Branch("v_pz", &m_v_pz);
  materialTree.Branch("v_phi", &m_v_phi);
  materialTree.Branch("v_eta", &m_v_eta);
  materialTree.Branch("t_X0", &m_tX0);
  materialTree.Branch("t_L0", &m_tL0);
  materialTree.Branch("mat_x", &m_step_x);
  materialTree.Branch("mat_y", &m_step_y);
  materialTree.Branch("mat_z", &m_step_z);
  materialTree.Branch("mat_r", &m_step_r);
  materialTree.Branch("mat_dx", &m_step_dx);
  materialTree.Branch("mat_dy", &m_step_dy);
  materialTree.Branch("mat_dz", &m_step_dz);
  materialTree.Branch("mat_step_length", &m_step_length);
  materialTree.Branch("mat_X0", &m_step_X0);
  materialTree.Branch("mat_L0", &m_step_L0);
  materialTree.Branch("mat_A", &m_step_A);
  materialTree.Branch("mat_Z", &m_step_Z);
  materialTree.Branch("mat_rho", &m_step_rho);

  if (m_prePostStepInfo) {
    materialTree.Branch("mat_sx", &m_step_sx);
    materialTree.Branch("mat_sy", &m_step_sy);
    materialTree.Branch("mat_sz", &m_step_sz);
    materialTree.Branch("mat_ex", &m_step_ex);
    materialTree.Branch("mat_ey", &m_step_ey);
    materialTree.Branch("mat_ez", &m_step_ez);
  }
  if (m_surfaceInfo) {
    materialTree.Branch("sur_id", &m_sur_id);
    materialTree.Branch("sur_x", &m_sur_x);
    materialTree.Branch("sur_y", &m_sur_y);
    materialTree.Branch("sur_z", &m_sur_z);
    materialTree.Branch("sur_r", &m_sur_r);
    materialTree.Branch("sur_distance", &m_sur_distance);
    materialTree.Branch("sur_pathCorrection", &m_sur_pathCorrection);
  }
  if (m_volumeInfo) {
    materialTree.Branch("vol_id", &m_vol_id);
  }
}

void Acts::RootMaterialTrackPayload::write(
    const GeometryContext& gctx, std::uint32_t eventNum,
    const RecordedMaterialTrack& materialTrack) {
  m_eventId = eventNum;
  // Clearing the vector first
  m_step_sx->clear();
  m_step_sy->clear();
  m_step_sz->clear();
  m_step_x->clear();
  m_step_y->clear();
  m_step_z->clear();
  m_step_r->clear();
  m_step_ex->clear();
  m_step_ey->clear();
  m_step_ez->clear();
  m_step_dx->clear();
  m_step_dy->clear();
  m_step_dz->clear();
  m_step_length->clear();
  m_step_X0->clear();
  m_step_L0->clear();
  m_step_A->clear();
  m_step_Z->clear();
  m_step_rho->clear();

  m_sur_id->clear();
  m_sur_x->clear();
  m_sur_y->clear();
  m_sur_z->clear();
  m_sur_r->clear();
  m_sur_distance->clear();
  m_sur_pathCorrection->clear();

  m_vol_id->clear();

  auto materialInteractions = materialTrack.second.materialInteractions;
  if (m_collapseInteractions) {
    std::vector<Acts::MaterialInteraction> collapsed;

    Acts::Vector3 positionSum = Acts::Vector3::Zero();
    double pathCorrectionSum = 0;

    for (std::size_t start = 0, end = 0; end < materialInteractions.size();
         ++end) {
      const auto& mintStart = materialInteractions[start];
      const auto& mintEnd = materialInteractions[end];

      positionSum += mintEnd.position;
      pathCorrectionSum += mintEnd.pathCorrection;

      const bool same =
          mintStart.materialSlab.material() == mintEnd.materialSlab.material();
      const bool last = end == materialInteractions.size() - 1;

      if (!same || last) {
        auto mint = mintStart;
        mint.position = positionSum / (end - start);
        mint.pathCorrection = pathCorrectionSum;

        collapsed.push_back(mint);

        start = end;
        positionSum = Acts::Vector3::Zero();
        pathCorrectionSum = 0;
      }
    }
    materialInteractions = std::move(collapsed);
  }

  // Reserve the vector then
  std::size_t mints = materialInteractions.size();
  m_step_sx->reserve(mints);
  m_step_sy->reserve(mints);
  m_step_sz->reserve(mints);
  m_step_x->reserve(mints);
  m_step_y->reserve(mints);
  m_step_z->reserve(mints);
  m_step_r->reserve(mints);
  m_step_ex->reserve(mints);
  m_step_ey->reserve(mints);
  m_step_ez->reserve(mints);
  m_step_dx->reserve(mints);
  m_step_dy->reserve(mints);
  m_step_dz->reserve(mints);
  m_step_length->reserve(mints);
  m_step_X0->reserve(mints);
  m_step_L0->reserve(mints);
  m_step_A->reserve(mints);
  m_step_Z->reserve(mints);
  m_step_rho->reserve(mints);

  m_sur_id->reserve(mints);
  m_sur_x->reserve(mints);
  m_sur_y->reserve(mints);
  m_sur_z->reserve(mints);
  m_sur_r->reserve(mints);
  m_sur_distance->reserve(mints);
  m_sur_pathCorrection->reserve(mints);

  m_vol_id->reserve(mints);

  // reset the global counter
  if (m_recalculateTotals) {
    m_tX0 = 0.;
    m_tL0 = 0.;
  } else {
    m_tX0 = materialTrack.second.materialInX0;
    m_tL0 = materialTrack.second.materialInL0;
  }

  // set the track information at vertex
  m_v_x = materialTrack.first.first.x();
  m_v_y = materialTrack.first.first.y();
  m_v_z = materialTrack.first.first.z();
  m_v_px = materialTrack.first.second.x();
  m_v_py = materialTrack.first.second.y();
  m_v_pz = materialTrack.first.second.z();
  m_v_phi = VectorHelpers::phi(materialTrack.first.second);
  m_v_eta = VectorHelpers::eta(materialTrack.first.second);

  // and now loop over the material
  for (const auto& mint : materialInteractions) {
    auto direction = mint.direction.normalized();

    // The material step position information
    m_step_x->push_back(mint.position.x());
    m_step_y->push_back(mint.position.y());
    m_step_z->push_back(mint.position.z());
    m_step_r->push_back(VectorHelpers::perp(mint.position));
    m_step_dx->push_back(direction.x());
    m_step_dy->push_back(direction.y());
    m_step_dz->push_back(direction.z());

    if (m_prePostStepInfo) {
      Acts::Vector3 prePos =
          mint.position - 0.5 * mint.pathCorrection * direction;
      Acts::Vector3 posPos =
          mint.position + 0.5 * mint.pathCorrection * direction;

      m_step_sx->push_back(prePos.x());
      m_step_sy->push_back(prePos.y());
      m_step_sz->push_back(prePos.z());
      m_step_ex->push_back(posPos.x());
      m_step_ey->push_back(posPos.y());
      m_step_ez->push_back(posPos.z());
    }

    // Store surface information
    if (m_surfaceInfo) {
      const Acts::Surface* surface = mint.surface;
      if (mint.intersectionID.value() != 0) {
        m_sur_id->push_back(mint.intersectionID.value());
        m_sur_pathCorrection->push_back(mint.pathCorrection);
        m_sur_x->push_back(mint.intersection.x());
        m_sur_y->push_back(mint.intersection.y());
        m_sur_z->push_back(mint.intersection.z());
        m_sur_r->push_back(VectorHelpers::perp(mint.intersection));
        m_sur_distance->push_back((mint.position - mint.intersection).norm());
      } else if (surface != nullptr) {
        auto sfIntersection =
            surface
                ->intersect(gctx, mint.position, mint.direction,
                            Acts::BoundaryTolerance::None())
                .closest();
        m_sur_id->push_back(surface->geometryId().value());
        m_sur_pathCorrection->push_back(1.0);
        m_sur_x->push_back(sfIntersection.position().x());
        m_sur_y->push_back(sfIntersection.position().y());
        m_sur_z->push_back(sfIntersection.position().z());
      } else {
        m_sur_id->push_back(Acts::GeometryIdentifier().value());
        m_sur_x->push_back(0);
        m_sur_y->push_back(0);
        m_sur_z->push_back(0);
        m_sur_pathCorrection->push_back(1.0);
      }
    }

    // store volume information
    if (m_volumeInfo) {
      Acts::GeometryIdentifier vlayerID;
      if (!mint.volume.empty()) {
        vlayerID = mint.volume.geometryId();
        m_vol_id->push_back(vlayerID.value());
      } else {
        vlayerID = vlayerID.withVolume(0)
                       .withBoundary(0)
                       .withLayer(0)
                       .withApproach(0)
                       .withSensitive(0);
        m_vol_id->push_back(vlayerID.value());
      }
    }

    // the material information
    const auto& mprops = mint.materialSlab;
    m_step_length->push_back(mprops.thickness());
    m_step_X0->push_back(mprops.material().X0());
    m_step_L0->push_back(mprops.material().L0());
    m_step_A->push_back(mprops.material().Ar());
    m_step_Z->push_back(mprops.material().Z());
    m_step_rho->push_back(mprops.material().massDensity());
    // re-calculate if defined to do so
    if (m_recalculateTotals) {
      m_tX0 += mprops.thicknessInX0();
      m_tL0 += mprops.thicknessInL0();
    }
  }
}

Acts::RecordedMaterialTrack Acts::RootMaterialTrackPayload::read() const {
  Acts::RecordedMaterialTrack rmTrack;
  // Fill the position and momentum
  rmTrack.first.first = Acts::Vector3(m_v_x, m_v_y, m_v_z);
  rmTrack.first.second = Acts::Vector3(m_v_px, m_v_py, m_v_pz);

  // Fill the individual steps
  std::size_t msteps = m_step_length->size();
  rmTrack.second.materialInteractions.reserve(msteps);
  rmTrack.second.materialInX0 = 0.;
  rmTrack.second.materialInL0 = 0.;

  for (std::size_t is = 0; is < msteps; ++is) {
    double s = (*m_step_length)[is];
    if (s == 0) {
      continue;
    }

    double mX0 = (*m_step_X0)[is];
    double mL0 = (*m_step_L0)[is];

    rmTrack.second.materialInX0 += s / mX0;
    rmTrack.second.materialInL0 += s / mL0;
    /// Fill the position & the material
    Acts::MaterialInteraction mInteraction;
    mInteraction.position =
        Acts::Vector3((*m_step_x)[is], (*m_step_y)[is], (*m_step_z)[is]);
    mInteraction.direction =
        Acts::Vector3((*m_step_dx)[is], (*m_step_dy)[is], (*m_step_dz)[is]);
    mInteraction.materialSlab = Acts::MaterialSlab(
        Acts::Material::fromMassDensity(mX0, mL0, (*m_step_A)[is],
                                        (*m_step_Z)[is], (*m_step_rho)[is]),
        s);
    if (m_surfaceInfo) {
      // add the surface information to the interaction this allows the
      // mapping to be speed up
      mInteraction.intersectionID = Acts::GeometryIdentifier((*m_sur_id)[is]);
      mInteraction.intersection =
          Acts::Vector3((*m_sur_x)[is], (*m_sur_y)[is], (*m_sur_z)[is]);
      mInteraction.pathCorrection = (*m_sur_pathCorrection)[is];
    } else {
      mInteraction.intersectionID = Acts::GeometryIdentifier();
      mInteraction.intersection = Acts::Vector3(0, 0, 0);
    }
    rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
  }
  return rmTrack;
}
