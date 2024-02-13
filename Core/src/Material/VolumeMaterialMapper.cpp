// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/VolumeMaterialMapper.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/VolumeCollector.hpp"
#include "Acts/Utilities/BinAdjustmentVolume.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Result.hpp"

#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace Acts {
struct EndOfWorldReached;
}  // namespace Acts

Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    const Config& cfg, StraightLineTGPropagator& propagator,
    std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_tgPropagator(
          std::make_shared<StraightLineTGPropagator>(std::move(propagator))),
      m_logger(std::move(slogger)) {}

Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    const Config& cfg, StraightLineDetPropagator& propagator,
    std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_detPropagator(
          std::make_shared<StraightLineDetPropagator>(std::move(propagator))),
      m_logger(std::move(slogger)) {}

std::unique_ptr<Acts::MaterialMappingState>
Acts::VolumeMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const TrackingGeometry& tGeometry) const {
  auto world = tGeometry.highestTrackingVolume();
  // The Volume material mapping state
  auto mState = std::make_unique<State>(gctx, mctx);
  resolveMaterialVolume(*mState, *world);
  return mState;
}

std::unique_ptr<Acts::MaterialMappingState>
Acts::VolumeMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const Experimental::Detector& detector) const {
  auto mState = std::make_unique<State>(gctx, mctx);
  for (const auto* volume : detector.volumes()) {
    resolveMaterialVolume(*mState, *volume);
  }
  return mState;
}

void Acts::VolumeMaterialMapper::resolveMaterialVolume(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for proto volume material.")

  if (tVolume.volumeMaterial() != nullptr) {
    checkAndInsert(mState, *tVolume.volumeMaterial(), tVolume.volumeBounds(),
                   tVolume.transform(), tVolume.geometryId());
  }

  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    ACTS_VERBOSE("- Check children volume ...");
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
  if (!tVolume.denseVolumes().empty()) {
    for (auto& sVolume : tVolume.denseVolumes()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::resolveMaterialVolume(
    State& mState, const Experimental::DetectorVolume& dVolume) const {
  ACTS_VERBOSE("Checking volume '" << dVolume.name()
                                   << "' for proto volume material.")

  if (dVolume.volumeMaterial() != nullptr) {
    checkAndInsert(mState, *dVolume.volumeMaterial(), dVolume.volumeBounds(),
                   dVolume.transform(mState.geoContext), dVolume.geometryId());
  }
  // Step down into the sub volume
  if (!dVolume.volumes().empty()) {
    ACTS_VERBOSE("- Check children volume ...");
    for (const auto* volume : dVolume.volumes()) {
      // Recursive call
      resolveMaterialVolume(mState, *volume);
    }
  }
}

void Acts::VolumeMaterialMapper::checkAndInsert(
    State& mState, const IVolumeMaterial& volumeMaterial,
    const VolumeBounds& volumeBounds, const Transform3& transform,
    const GeometryIdentifier& geoID) const {
  std::size_t volumeID = geoID.volume();
  ACTS_DEBUG("Material volume found with volumeID " << volumeID);
  ACTS_DEBUG("       - ID is " << geoID);

  auto psm = dynamic_cast<const ProtoVolumeMaterial*>(&volumeMaterial);
  // Get the bin utility: try proxy material first
  const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
  if (bu != nullptr) {
    // Screen output for Binned Proto material
    ACTS_DEBUG("       - (proto) binning is " << *bu);
    // Now update
    BinUtility buAdjusted = adjustBinUtility(*bu, volumeBounds, transform);
    // Screen output for Binned Proto material
    ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
    mState.materialBin[geoID] = buAdjusted;
    if (bu->dimensions() == 0) {
      ACTS_DEBUG("Binning of dimension 0 create AccumulatedVolumeMaterial");
      Acts::AccumulatedVolumeMaterial homogeneousAccumulation;
      mState.homogeneousGrid[geoID] = homogeneousAccumulation;
    } else if (bu->dimensions() == 2) {
      ACTS_DEBUG("Binning of dimension 2 create 2D Grid");
      std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
      Acts::Grid2D Grid = createGrid2D(buAdjusted, transfoGlobalToLocal);
      mState.grid2D.insert(std::make_pair(geoID, Grid));
      mState.transform2D.insert(std::make_pair(geoID, transfoGlobalToLocal));
    } else if (bu->dimensions() == 3) {
      ACTS_DEBUG("Binning of dimension 3 create 3D Grid");
      std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
      Acts::Grid3D Grid = createGrid3D(buAdjusted, transfoGlobalToLocal);
      mState.grid3D.insert(std::make_pair(geoID, Grid));
      mState.transform3D.insert(std::make_pair(geoID, transfoGlobalToLocal));
    } else {
      throw std::invalid_argument(
          "Incorrect bin dimension, only 0, 2 and 3 are accepted");
    }
    return;
  }
  // Second attempt: 2D binned material
  auto bmp2 = dynamic_cast<
      const InterpolatedMaterialMap<MaterialMapper<Acts::MaterialGrid2D>>*>(
      &volumeMaterial);
  bu = (bmp2 != nullptr) ? (&bmp2->binUtility()) : nullptr;
  if (bu != nullptr) {
    // Screen output for Binned Proto material
    ACTS_DEBUG("       - (2D grid) binning is " << *bu);
    mState.materialBin[geoID] = *bu;
    std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid2D Grid = createGrid2D(*bu, transfoGlobalToLocal);
    mState.grid2D.insert(std::make_pair(geoID, Grid));
    mState.transform2D.insert(std::make_pair(geoID, transfoGlobalToLocal));
    return;
  }
  // Third attempt: 3D binned material
  auto bmp3 = dynamic_cast<
      const InterpolatedMaterialMap<MaterialMapper<Acts::MaterialGrid3D>>*>(
      &volumeMaterial);
  bu = (bmp3 != nullptr) ? (&bmp3->binUtility()) : nullptr;
  if (bu != nullptr) {
    // Screen output for Binned Proto material
    ACTS_DEBUG("       - (3D grid) binning is " << *bu);
    mState.materialBin[geoID] = *bu;
    std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid3D Grid = createGrid3D(*bu, transfoGlobalToLocal);
    mState.grid3D.insert(std::make_pair(geoID, Grid));
    mState.transform3D.insert(std::make_pair(geoID, transfoGlobalToLocal));
    return;
  } else {
    // Create a homogeneous type of material
    ACTS_DEBUG("       - this is homogeneous material.");
    BinUtility buHomogeneous;
    mState.materialBin[geoID] = buHomogeneous;
    Acts::AccumulatedVolumeMaterial homogeneousAccumulation;
    mState.homogeneousGrid[geoID] = homogeneousAccumulation;
    return;
  }
}

void Acts::VolumeMaterialMapper::createExtraHits(
    State& mState,
    std::pair<const GeometryIdentifier, BinUtility>& currentBinning,
    Acts::MaterialSlab properties, const Vector3& position,
    Vector3 direction) const {
  if (currentBinning.second.dimensions() == 0) {
    // Writing homogeneous material for the current volumes no need to create
    // extra hits. We directly accumulate the material
    mState.homogeneousGrid[currentBinning.first].accumulate(properties);
    return;
  }

  // Computing the extra hits properties based on the mappingStep length
  int volumeStep =
      static_cast<int>(std::floor(properties.thickness() / m_cfg.mappingStep));
  float remainder = properties.thickness() - m_cfg.mappingStep * volumeStep;
  properties.scaleThickness(m_cfg.mappingStep / properties.thickness());
  direction = direction * (m_cfg.mappingStep / direction.norm());

  for (int extraStep = 0; extraStep < volumeStep; extraStep++) {
    Vector3 extraPosition = position + extraStep * direction;
    // Create additional extrapolated points for the grid mapping

    if (currentBinning.second.dimensions() == 2) {
      auto grid = mState.grid2D.find(currentBinning.first);
      if (grid != mState.grid2D.end()) {
        // Find which grid bin the material fall into then accumulate
        Acts::Grid2D::index_t index = grid->second.localBinsFromLowerLeftEdge(
            mState.transform2D[currentBinning.first](extraPosition));
        grid->second.atLocalBins(index).accumulate(properties);
      } else {
        throw std::domain_error("No grid 2D was found");
      }
    } else if (currentBinning.second.dimensions() == 3) {
      auto grid = mState.grid3D.find(currentBinning.first);
      if (grid != mState.grid3D.end()) {
        // Find which grid bin the material fall into then accumulate
        Acts::Grid3D::index_t index = grid->second.localBinsFromLowerLeftEdge(
            mState.transform3D[currentBinning.first](extraPosition));
        grid->second.atLocalBins(index).accumulate(properties);
      } else {
        throw std::domain_error("No grid 3D was found");
      }
    }
  }

  if (remainder > 0) {
    // We need to had an additional extra hit with the remainder length. Adjust
    // the thickness of the last extrapolated step
    properties.scaleThickness(remainder / properties.thickness());
    Vector3 extraPosition = position + volumeStep * direction;
    if (currentBinning.second.dimensions() == 2) {
      auto grid = mState.grid2D.find(currentBinning.first);
      if (grid != mState.grid2D.end()) {
        // Find which grid bin the material fall into then accumulate
        Acts::Grid2D::index_t index = grid->second.localBinsFromLowerLeftEdge(
            mState.transform2D[currentBinning.first](extraPosition));
        grid->second.atLocalBins(index).accumulate(properties);
      } else {
        throw std::domain_error("No grid 2D was found");
      }
    } else if (currentBinning.second.dimensions() == 3) {
      auto grid = mState.grid3D.find(currentBinning.first);
      if (grid != mState.grid3D.end()) {
        // Find which grid bin the material fall into then accumulate
        Acts::Grid3D::index_t index = grid->second.localBinsFromLowerLeftEdge(
            mState.transform3D[currentBinning.first](extraPosition));
        grid->second.atLocalBins(index).accumulate(properties);
      } else {
        throw std::domain_error("No grid 3D was found");
      }
    }
  }
}

Acts::MaterialMappingResult Acts::VolumeMaterialMapper::finalizeMaps(
    MaterialMappingState& mState) const {
  State& state = static_cast<State&>(mState);

  // Fill the material maps
  Acts::MaterialMappingResult maps{};

  // iterate over the volumes
  for (auto& matBin : state.materialBin) {
    ACTS_DEBUG("Create the material for volume  " << matBin.first);
    if (matBin.second.dimensions() == 0) {
      // Average the homogeneous volume material then store it
      ACTS_DEBUG("Homogeneous material volume");
      Acts::Material mat = state.homogeneousGrid[matBin.first].average();
      maps.volumeMaterial[matBin.first] =
          std::make_unique<HomogeneousVolumeMaterial>(mat);
    } else if (matBin.second.dimensions() == 2) {
      // Average the material in the 2D grid then create an
      // InterpolatedMaterialMap
      ACTS_DEBUG("Grid material volume");
      auto grid = state.grid2D.find(matBin.first);
      if (grid != state.grid2D.end()) {
        Acts::MaterialGrid2D matGrid = mapMaterialPoints(grid->second);
        MaterialMapper<Acts::MaterialGrid2D> matMap(
            state.transform2D[matBin.first], matGrid);
        maps.volumeMaterial[matBin.first] = std::make_unique<
            InterpolatedMaterialMap<MaterialMapper<Acts::MaterialGrid2D>>>(
            std::move(matMap), matBin.second);
      } else {
        throw std::domain_error("No grid 2D was found");
      }
    } else if (matBin.second.dimensions() == 3) {
      // Average the material in the 3D grid then create an
      // InterpolatedMaterialMap
      ACTS_DEBUG("Grid material volume");
      auto grid = state.grid3D.find(matBin.first);
      if (grid != state.grid3D.end()) {
        Acts::MaterialGrid3D matGrid = mapMaterialPoints(grid->second);
        MaterialMapper<Acts::MaterialGrid3D> matMap(
            state.transform3D[matBin.first], matGrid);
        maps.volumeMaterial[matBin.first] = std::make_unique<
            InterpolatedMaterialMap<MaterialMapper<Acts::MaterialGrid3D>>>(
            std::move(matMap), matBin.second);
      } else {
        throw std::domain_error("No grid 3D was found");
      }
    } else {
      throw std::invalid_argument(
          "Incorrect bin dimension, only 0, 2 and 3 are accepted");
    }
  }
  return maps;
}

std::array<Acts::RecordedMaterialTrack, 2u>
Acts::VolumeMaterialMapper::mapMaterialTrack(
    MaterialMappingState& mState, const RecordedMaterialTrack& mTrack) const {
  using VectorHelpers::makeVector4;

  State& state = static_cast<State&>(mState);

  RecordedMaterialTrack mapped{mTrack.first, RecordedMaterial{}};
  RecordedMaterialTrack unmapped{mTrack.first, RecordedMaterial{}};

  // Neutral curvilinear parameters
  NeutralCurvilinearTrackParameters start(
      makeVector4(mTrack.first.first, 0), mTrack.first.second,
      1 / mTrack.first.second.norm(), std::nullopt,
      NeutralParticleHypothesis::geantino());

  // Prepare Action list and abort list
  using MaterialVolumeCollector = VolumeCollector<MaterialVolumeSelector>;
  using ActionList = ActionList<MaterialVolumeCollector>;
  using AbortList = AbortList<EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options(state.geoContext,
                                                   state.magFieldContext);

  // Now collect the material layers by using the straight line propagator
  MaterialVolumeCollector::result_type mcResult;
  // Now collect the material layers by using the straight line propagator
  if (m_tgPropagator) {
    const auto& result = m_tgPropagator->propagate(start, options).value();
    mcResult = result.get<MaterialVolumeCollector::result_type>();
  } else {
    // Now collect the material layers by using the straight line propagator
    const auto& result = m_detPropagator->propagate(start, options).value();
    mcResult = result.get<MaterialVolumeCollector::result_type>();
  }
  auto mappingVolumes = mcResult.collected;

  // Retrieve the recorded material from the recorded material track
  auto& rMaterial = mTrack.second.materialInteractions;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material steps to map.")

  // These should be mapped onto the volumes
  ACTS_VERBOSE("Found     " << mappingVolumes.size()
                            << " mapping volumes for this track.");
  ACTS_VERBOSE("Mapping volumes are :")
  for (auto& mVolumes : mappingVolumes) {
    ACTS_VERBOSE(" - Volume : " << mVolumes.geometryId() << " at position = ("
                                << mVolumes.position.x() << ", "
                                << mVolumes.position.y() << ", "
                                << mVolumes.position.z() << ")");

    // mappingVolumes.push_back(mVolumes);
  }
  // Run the mapping process, i.e. take the recorded material and map it
  // onto the mapping volume:
  auto rmIter = rMaterial.begin();
  auto volIter = mappingVolumes.begin();

  // Use those to minimize the lookup
  GeometryIdentifier lastID = GeometryIdentifier();
  GeometryIdentifier currentID = GeometryIdentifier();
  auto currentBinning = state.materialBin.end();

  // store end position of the last material slab
  Acts::Vector3 lastPositionEnd = {0, 0, 0};
  Acts::Vector3 direction = {0, 0, 0};

  if (volIter != mappingVolumes.end()) {
    lastPositionEnd = volIter->position;
  }

  // loop over all the material hit in the track or until there no more volume
  // to map onto
  while (rmIter != rMaterial.end() && volIter != mappingVolumes.end()) {
    // Check if the the recorded material interaction is vetoed
    if (m_cfg.veto(*rmIter)) {
      unmapped.second.materialInX0 += rmIter->materialSlab.thicknessInX0();
      unmapped.second.materialInL0 += rmIter->materialSlab.thicknessInL0();
      unmapped.second.materialInteractions.push_back(*rmIter);
      ++rmIter;
      continue;
    }
    // Assume it's mapped then
    auto mappedRecord = (*rmIter);
    if (volIter != mappingVolumes.end() &&
        !volIter->inside(state.geoContext, rmIter->position)) {
      // Check if the material point is past the entry point to the current
      // volume (this prevent switching volume before the first volume has been
      // reached)
      double distVol = (volIter->position - mTrack.first.first).norm();
      double distMat = (rmIter->position - mTrack.first.first).norm();
      if (distMat - distVol > s_epsilon) {
        // Switch to next material volume
        ++volIter;
        continue;
      }
    }
    if (volIter != mappingVolumes.end() &&
        volIter->inside(state.geoContext, rmIter->position, s_epsilon)) {
      currentID = volIter->geometryId();
      direction = rmIter->direction;
      if (!(currentID == lastID)) {
        // Let's (re-)assess the information
        lastID = currentID;
        lastPositionEnd = volIter->position;
        currentBinning = state.materialBin.find(currentID);
      }
      // If the current volume has a ProtoVolumeMaterial
      // and the material hit has a non 0 thickness
      if (currentBinning != state.materialBin.end() &&
          rmIter->materialSlab.thickness() > 0) {
        // check if there is vacuum between this material point and the last one
        float vacuumThickness = (rmIter->position - lastPositionEnd).norm();
        if (vacuumThickness > s_epsilon) {
          auto properties = Acts::MaterialSlab(vacuumThickness);
          // Creat vacuum hits
          createExtraHits(state, *currentBinning, properties, lastPositionEnd,
                          direction);
        }
        // Determine the position of the last material slab using the track
        // direction
        direction =
            direction * (rmIter->materialSlab.thickness() / direction.norm());
        lastPositionEnd = rmIter->position + direction;
        // create additional material point
        createExtraHits(state, *currentBinning, rmIter->materialSlab,
                        rmIter->position, direction);
      }
      // mappedRecord.volume = volIter->volume;
      mappedRecord.intersectionID = currentID;
      mappedRecord.intersection = rmIter->position;
      mapped.second.materialInX0 += rmIter->materialSlab.thicknessInX0();
      mapped.second.materialInL0 += rmIter->materialSlab.thicknessInL0();
      mapped.second.materialInteractions.push_back(*rmIter);
    } else {
      unmapped.second.materialInX0 += rmIter->materialSlab.thicknessInX0();
      unmapped.second.materialInL0 += rmIter->materialSlab.thicknessInL0();
      unmapped.second.materialInteractions.push_back(*rmIter);
      ++rmIter;
      continue;
    }
    ++rmIter;
  }

  return {mapped, unmapped};
}
