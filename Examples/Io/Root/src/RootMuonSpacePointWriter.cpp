// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMuonSpacePointWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::VectorHelpers;
/// @brief Converts a surface Identifier to the one of the surrounding volume
/// @param id: Surface identifier to convert
constexpr GeometryIdentifier toChamberId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume()).withLayer(id.layer());
}
/// @brief Pushes value back to a vector and applies a static cast at the same
/// @param vec: Vector into which the value is pushed back
/// @param val: Value to push
template <typename T, typename T1>
void castPush(std::vector<T>& vec, const T1& val) {
  // MARK: fpeMaskBegin(FLTUND, 1, #-1)
  const T castedVal = static_cast<T>(val);
  vec.push_back(castedVal);
  // MARK: fpeMaskEnd(FLTUND)
}
/// @brief Converts an angle in radians into an angle in degree
constexpr double inDeg(const double radians) {
  if (Acts::abs(radians) < std::numeric_limits<float>::epsilon()) {
    return 0.;
  }
  using namespace Acts::UnitLiterals;
  return radians / 1._degree;
}
namespace ActsExamples {
RootMuonSpacePointWriter::RootMuonSpacePointWriter(const Config& config,
                                                   Logging::Level level)
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

  m_tree->Branch("event_id", &m_eventId);
  m_tree->Branch("spacePoint_bucketId", &m_bucketId);
  m_tree->Branch("spacePoint_geometryId", &m_geometryId);
  m_tree->Branch("spacePoint_muonId", &m_muonId);
  m_tree->Branch("spacePoint_localPosX", &m_localPositionX);
  m_tree->Branch("spacePoint_localPosY", &m_localPositionY);
  m_tree->Branch("spacePoint_localPosZ", &m_localPositionZ);

  m_tree->Branch("spacePoint_sensorDirTheta", &m_sensorDirectionTheta);
  m_tree->Branch("spacePoint_sensorDirPhi", &m_sensorDirectionPhi);

  m_tree->Branch("spacePoint_toNextDirTheta", &m_toNextSensorTheta);
  m_tree->Branch("spacePoint_toNextDirPhi", &m_toNextSensorPhi);

  m_tree->Branch("spacePoint_covLoc0", &m_covLoc0);
  m_tree->Branch("spacePoint_covLoc1", &m_covLoc1);
  m_tree->Branch("spacePoint_covT", &m_covLocT);

  m_tree->Branch("spacePoint_driftRadius", &m_driftR);
  m_tree->Branch("spacePoint_time", &m_time);

  if (m_cfg.writeGlobal && m_cfg.trackingGeometry == nullptr) {
    throw std::runtime_error(
        "RootMuonSpacePointWriter() - Global coordinates can only be written "
        "with a tracking geometry");
  }
  if (m_cfg.writeGlobal) {
    m_tree->Branch("global_positionX", &m_globalPosX);
    m_tree->Branch("global_positionY", &m_globalPosY);
    m_tree->Branch("global_positionZ", &m_globalPosZ);

    m_tree->Branch("global_lowEdgeX", &m_lowEdgeX);
    m_tree->Branch("global_lowEdgeY", &m_lowEdgeY);
    m_tree->Branch("global_lowEdgeZ", &m_lowEdgeZ);

    m_tree->Branch("global_highEdgeX", &m_highEdgeX);
    m_tree->Branch("global_highEdgeY", &m_highEdgeY);
    m_tree->Branch("global_highEdgeZ", &m_highEdgeZ);
  }
}

RootMuonSpacePointWriter::~RootMuonSpacePointWriter() = default;
ProcessCode RootMuonSpacePointWriter::finalize() {
  m_file->cd();
  m_file->Write();
  m_file.reset();
  ACTS_INFO("Wrote muon spacepoints to tree '" << m_cfg.treeName << "' in '"
                                               << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}
ProcessCode RootMuonSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const MuonSpacePointContainer& hits) {
  std::lock_guard lock{m_mutex};
  m_eventId = ctx.eventNumber;
  const Acts::GeometryContext gctx{};
  for (const auto& [counter, bucket] : enumerate(hits)) {
    for (const MuonSpacePoint& writeMe : bucket) {
      ACTS_VERBOSE("Dump space point " << writeMe);

      castPush(m_bucketId, counter);
      castPush(m_geometryId, writeMe.geometryId().value());
      castPush(m_muonId, writeMe.id().toInt());
      castPush(m_localPositionX, writeMe.localPosition().x());
      castPush(m_localPositionY, writeMe.localPosition().y());
      castPush(m_localPositionZ, writeMe.localPosition().z());

      castPush(m_sensorDirectionTheta, inDeg(theta(writeMe.sensorDirection())));
      castPush(m_sensorDirectionPhi, inDeg(phi(writeMe.sensorDirection())));

      castPush(m_toNextSensorTheta, inDeg(theta(writeMe.toNextSensor())));
      castPush(m_toNextSensorPhi, inDeg(phi(writeMe.toNextSensor())));

      const auto& cov = writeMe.covariance();
      {
        using namespace Experimental::detail;
        using enum CompSpacePointAuxiliaries::ResidualIdx;
        castPush(m_covLoc0, cov[toUnderlying(nonBending)]);
        castPush(m_covLoc1, cov[toUnderlying(bending)]);
        castPush(m_covLocT, cov[toUnderlying(time)]);
      }
      castPush(m_driftR, writeMe.driftRadius());
      castPush(m_time, writeMe.time());
      if (!m_cfg.writeGlobal || writeMe.geometryId() == GeometryIdentifier{}) {
        continue;
      }
      const Surface* surface =
          m_cfg.trackingGeometry->findSurface(writeMe.geometryId());
      assert(surface != nullptr);
      const TrackingVolume* chambVol =
          m_cfg.trackingGeometry->findVolume(toChamberId(writeMe.geometryId()));
      assert(chambVol != nullptr);

      const Vector3 globPos = chambVol->transform() *
                              AngleAxis3{-90._degree, Vector3::UnitZ()} *
                              writeMe.localPosition();
      castPush(m_globalPosX, globPos.x());
      castPush(m_globalPosY, globPos.y());
      castPush(m_globalPosZ, globPos.z());

      const auto& bounds{surface->bounds()};
      const auto& trf{surface->localToGlobal(gctx)};
      Acts::Vector3 lowEdge{Vector3::Zero()};
      Acts::Vector3 highEdge{Vector3::Zero()};
      switch (bounds.type()) {
        using enum Acts::SurfaceBounds::BoundsType;
        case eLine: {
          const auto& lBounds = static_cast<const LineBounds&>(bounds);
          const double l = lBounds.get(LineBounds::eHalfLengthZ);
          lowEdge = trf * (-l * Vector3::UnitZ());
          highEdge = trf * (l * Vector3::UnitZ());
          break;
        }
        case eRectangle: {
          const auto& rBounds = static_cast<const RectangleBounds&>(bounds);
          const double l =
              rBounds.get(writeMe.measuresLoc1() ? RectangleBounds::eMaxX
                                                 : RectangleBounds::eMaxY);
          const auto dimIdx =
              static_cast<std::uint32_t>(writeMe.measuresLoc1() == false);
          lowEdge = trf * (-l * Vector3::Unit(dimIdx));
          highEdge = trf * (l * Vector3::Unit(dimIdx));
          break;
        }
        case eTrapezoid: {
          ACTS_WARNING(__FILE__ << ":" << __LINE__ << " Implement me");
          break;
        }
        /// Other surface bound types are not supported
        default:
          ACTS_ERROR("Unsupported bounds " << surface->toString(gctx));
          return ProcessCode::ABORT;
      }
      castPush(m_lowEdgeX, lowEdge.x());
      castPush(m_lowEdgeY, lowEdge.y());
      castPush(m_lowEdgeZ, lowEdge.z());
      castPush(m_highEdgeX, highEdge.x());
      castPush(m_highEdgeY, highEdge.y());
      castPush(m_highEdgeZ, highEdge.z());
    }
  }
  m_tree->Fill();

  m_geometryId.clear();
  m_bucketId.clear();
  m_muonId.clear();
  m_localPositionX.clear();
  m_localPositionY.clear();
  m_localPositionZ.clear();
  m_sensorDirectionTheta.clear();
  m_sensorDirectionPhi.clear();
  m_toNextSensorTheta.clear();
  m_toNextSensorPhi.clear();
  m_covLoc0.clear();
  m_covLoc1.clear();
  m_covLocT.clear();
  m_driftR.clear();
  m_time.clear();

  m_globalPosX.clear();
  m_globalPosY.clear();
  m_globalPosZ.clear();

  m_lowEdgeX.clear();
  m_lowEdgeY.clear();
  m_lowEdgeZ.clear();

  m_highEdgeX.clear();
  m_highEdgeY.clear();
  m_highEdgeZ.clear();
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
