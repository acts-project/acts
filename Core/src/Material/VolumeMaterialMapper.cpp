// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Material/MaterialMapUtils.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/BinAdjustmentVolume.hpp"

namespace {
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid2D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2) {
  // get the number of bins
  size_t nBinsAxis1 = gridAxis1[2];
  size_t nBinsAxis2 = gridAxis2[2];

  // get the minimum and maximum
  double minAxis1 = gridAxis1[0];
  double minAxis2 = gridAxis2[0];
  double maxAxis1 = gridAxis1[1];
  double maxAxis2 = gridAxis2[1];
  // calculate maxima (add one last bin, because bin value always corresponds
  // to
  // left boundary)
  double stepAxis1 = std::fabs(maxAxis1 - minAxis1) / (nBinsAxis1 - 1);
  double stepAxis2 = std::fabs(maxAxis2 - minAxis2) / (nBinsAxis2 - 1);
  maxAxis1 += stepAxis1;
  maxAxis2 += stepAxis2;

  // Create the axis for the grid
  EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);

  // The material mapping grid
  return Grid2D(std::make_tuple(std::move(axis1), std::move(axis2)));
}

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @param [in] gridAxis3 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid3D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2,
                  std::array<double, 3> gridAxis3) {
  // get the number of bins
  size_t nBinsAxis1 = gridAxis1[2];
  size_t nBinsAxis2 = gridAxis2[2];
  size_t nBinsAxis3 = gridAxis3[2];

  // get the minimum and maximum
  double minAxis1 = gridAxis1[0];
  double minAxis2 = gridAxis2[0];
  double minAxis3 = gridAxis3[0];
  double maxAxis1 = gridAxis1[1];
  double maxAxis2 = gridAxis2[1];
  double maxAxis3 = gridAxis3[1];
  // calculate maxima (add one last bin, because bin value always corresponds
  // to
  // left boundary)
  double stepAxis1 = std::fabs(maxAxis1 - minAxis1) / (nBinsAxis1 - 1);
  double stepAxis2 = std::fabs(maxAxis2 - minAxis2) / (nBinsAxis2 - 1);
  double stepAxis3 = std::fabs(maxAxis3 - minAxis3) / (nBinsAxis3 - 1);
  maxAxis1 += stepAxis1;
  maxAxis2 += stepAxis2;
  maxAxis3 += stepAxis3;

  // Create the axis for the grid
  EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);
  EAxis axis3(minAxis3, maxAxis3, nBinsAxis3);

  // The material mapping grid
  return Grid3D(
      std::make_tuple(std::move(axis1), std::move(axis2), std::move(axis3)));
}

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const Acts::RecordedMaterialPoint& mPoints,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid2D::index_t index = matchToGridPoint(rm.second, grid);
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid2D::point_t min = grid.minPosition();
  Grid2D::point_t max = grid.maxPosition();
  Grid2D::index_t nBins = grid.numLocalBins();

  EAxis axis1(min[0], max[0], nBins[0]);
  EAxis axis2(min[1], max[1], nBins[1]);

  // Build the grid and fill it with data
  MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().classificationNumbers();
  }

  return mGrid;
}

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid3D mapMaterialPoints(
    Grid3D& grid, const Acts::RecordedMaterialPoint& mPoints,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid3D::index_t index = matchToGridPoint(rm.second, grid);
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid3D::point_t min = grid.minPosition();
  Grid3D::point_t max = grid.maxPosition();
  Grid3D::index_t nBins = grid.numLocalBins();

  EAxis axis1(min[0], max[0], nBins[0]);
  EAxis axis2(min[1], max[1], nBins[1]);
  EAxis axis3(min[2], max[2], nBins[2]);

  // Build the grid and fill it with data
  MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().classificationNumbers();
  }
  return mGrid;
}
}  // namespace

MaterialGrid2D Acts::createMaterialGrid(
    std::array<double, 3> gridAxis1, std::array<double, 3> gridAxis2,
    const RecordedMaterialPoint& mPoints,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint) {
  Grid2D grid = createGrid(std::move(gridAxis1), std::move(gridAxis2));
  return mapMaterialPoints(grid, mPoints, matchToGridPoint);
}

MaterialGrid3D Acts::createMaterialGrid(
    std::array<double, 3> gridAxis1, std::array<double, 3> gridAxis2,
    std::array<double, 3> gridAxis3, const RecordedMaterialPoint& mPoints,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint) {
  Grid3D grid = createGrid(std::move(gridAxis1), std::move(gridAxis2),
                           std::move(gridAxis3));
  return mapMaterialPoints(grid, mPoints, matchToGridPoint);
}

/// @brief Rough searcher for closest point
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid3D::index_t mapMaterial3D(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0, indexZ = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      for (size_t k = 1; k < grid.numLocalBins()[2]+1; k++) {
        // Search the closest distance - elements are ordered
        double dX = grid.lowerLeftBinEdge({{i, j, k}})[0] - matPos.x();
        double dY = grid.lowerLeftBinEdge({{i, j, k}})[1] - matPos.y();
        double dZ = grid.lowerLeftBinEdge({{i, j, k}})[2] - matPos.z();

        if (std::sqrt(dX * dX + dY * dY + dZ * dZ) < dist) {
          // Store distance and index
          dist = std::sqrt(dX * dX + dY * dY + dZ * dZ);
          indexX = i;
          indexY = j;
          indexZ = k;
        } else {  // Break if distance becomes larger
          break;
        }
      }
    }
  }
  return {{indexX, indexY, indexZ}};
}

/// @brief Rough searcher for closest point in cylindrical coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid3D::index_t mapMaterialCylinder(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0, indexZ = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      for (size_t k = 1; k < grid.numLocalBins()[2]+1; k++) {
        double X = grid.lowerLeftBinEdge({{i, j, k}})[0]*cos(grid.lowerLeftBinEdge({{i, j, k}})[1]);
        double Y = grid.lowerLeftBinEdge({{i, j, k}})[0]*sin(grid.lowerLeftBinEdge({{i, j, k}})[1]);
        // Search the closest distance - elements are ordered
        double dX = X - matPos.x();
        double dY = Y - matPos.y();
        double dZ = grid.lowerLeftBinEdge({{i, j, k}})[2] - matPos.z();
        if (std::sqrt(dX * dX + dY * dY + dZ * dZ) < dist) {
          // Store distance and index
          dist = std::sqrt(dX * dX + dY * dY + dZ * dZ);
          indexX = i;
          indexY = j;
          indexZ = k;
        } else {  // Break if distance becomes larger
          // break;
        }
      }
    }
  }
  return {{indexX, indexY, indexZ}};
}

/// @brief Rough searcher for closest point
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
void updateMapOverfow(MaterialGrid3D& grid, const Acts::BinUtility& bu){
  auto& binningData = bu.binningData();
  size_t index1 = 0, index2 = 0, index3 = 0;
  size_t lastBin = 0;
  Grid3D::index_t localOver = {{0,0,0}};
  Grid3D::index_t localCopy = {{0,0,0}};
  for (index1 = 0; index1 < 3; index1++) {
    if(index1==0){
      index2 = 1;
      index3 = 2;
    }
    else if(index1==1){
      index2 = 0;
      index3 = 2;
    }
    else if(index1==2){
      index2 = 0;
      index3 = 1;
    }
    lastBin = grid.numLocalBins()[index1]+1;
    for (size_t j = 0; j < grid.numLocalBins()[index2]+1; j++) {
      for (size_t k = 0; k < grid.numLocalBins()[index3]+1; k++) {
        // Update the Underflow of the map
        localOver.at(index1) = 0; localOver.at(index2) = j; localOver.at(index3) = k;
        if(binningData[index1].option == Acts::closed){
          localCopy.at(index1) = lastBin-1; localCopy.at(index2) = j; localCopy.at(index3) = k;
          grid.atLocalBins(localOver) = grid.atLocalBins(localCopy);
        }
        else{
          localCopy.at(index1) = 1; localCopy.at(index2) = j; localCopy.at(index3) = k;
          grid.atLocalBins(localOver) = grid.atLocalBins(localCopy);
        }
        // Update the Overflow of the map
        localOver.at(index1) = lastBin; localOver.at(index2) = j; localOver.at(index3) = k;
        if(binningData[index1].option == Acts::closed){
          localCopy.at(index1) = 1; localCopy.at(index2) = j; localCopy.at(index3) = k;
          grid.atLocalBins(localOver) = grid.atLocalBins(localCopy);
        }
        else{
          localCopy.at(index1) = lastBin-1; localCopy.at(index2) = j; localCopy.at(index3) = k;
          grid.atLocalBins(localOver) = grid.atLocalBins(localCopy);
        }
      }
    }
  }
}


Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    const Config& cfg, StraightLinePropagator propagator,
     std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_propagator(std::move(propagator)),
      m_logger(std::move(slogger)) {}

Acts::VolumeMaterialMapper::State Acts::VolumeMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const TrackingGeometry& tGeometry) const {
  // Parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();

  // The Surface material mapping state
  State mState(gctx, mctx);
  resolveMaterialVolume(mState, *world);
  collectMaterialSurface(mState, *world);

  return mState;
}

void Acts::VolumeMaterialMapper::resolveMaterialVolume(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- Insert Volume ...");
  checkAndInsert(mState, tVolume);

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

void Acts::VolumeMaterialMapper::checkAndInsert(State& mState,
                                                 const TrackingVolume& volume) const
                                                 {
  auto volumeMaterial = volume.volumeMaterial();
  // Check if the volume has a proxy
  if (volumeMaterial != nullptr) {

    auto geoID = volume.geoID();
    size_t volumeID = geoID.volume();
    ACTS_DEBUG("Material volume found with volumeID " << volumeID);
    ACTS_DEBUG("       - ID is " << geoID);

    RecordedMaterialPoint mat;
    mState.recordedMaterial[geoID] = mat;

    // We need a dynamic_cast to either a volume material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoVolumeMaterial*>(volumeMaterial);
    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, volume);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
      mState.materialBin[geoID] = buAdjusted;
      return;
    }
    //TODO!!!
    // Second attempt: binned material
    auto bmp = dynamic_cast<const InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>*>(volumeMaterial);
    bu = (bmp != nullptr) ? (&bmp->binUtility()) : nullptr;
    // Creaete a binned type of material
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning is " << *bu);
      mState.materialBin[geoID] = *bu;
      return;
    } else {
      // Create a homogeneous type of material
      ACTS_DEBUG("       - this is homogeneous material.");
      return;
    }
  }
}

void Acts::VolumeMaterialMapper::collectMaterialSurface(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- boundary surfaces ...");
  // Check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces()) {
    if(bSurface->surfaceRepresentation().surfaceMaterial() != nullptr){
      mState.surfaceMaterial[bSurface->surfaceRepresentation().geoID()] = bSurface->surfaceRepresentation().surfaceMaterialSharedPtr();
    }
  }

  ACTS_VERBOSE("- confined layers ...");
  // Check the confined layers
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // Take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // Check the representing surface
        if(cLayer->surfaceRepresentation().surfaceMaterial() != nullptr){
          mState.surfaceMaterial[cLayer->surfaceRepresentation().geoID()] = cLayer->surfaceRepresentation().surfaceMaterialSharedPtr();
        }
        // Get the approach surfaces if present
        if (cLayer->approachDescriptor() != nullptr) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces()) {
            if (aSurface != nullptr) {
              if(aSurface->surfaceMaterial() != nullptr){
                mState.surfaceMaterial[aSurface->geoID()] = aSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
        // Get the sensitive surface is present
        if (cLayer->surfaceArray() != nullptr) {
          // Sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->surfaces()) {
            if (sSurface != nullptr) {
              if(sSurface->surfaceMaterial() != nullptr){
                mState.surfaceMaterial[sSurface->geoID()] = sSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
      }
    }
  }
  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      collectMaterialSurface(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::finalizeMaps(State& mState) const {
  // iterate over the map to call the total average
  for (auto& recMaterial : mState.recordedMaterial) {
    ACTS_DEBUG("Create the material for volume  " << recMaterial.first);
    auto bu = mState.materialBin[recMaterial.first].binningData();
    std::array<double, 3> gridAxis1;
    std::array<double, 3> gridAxis2;
    std::array<double, 3> gridAxis3;
    bool iscylinder = false;
    bool iscube = false;
    std::function<Acts::Vector3D(Acts::Vector3D)> transfoGlobalToLocal;
    std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)> mapMaterial;
    for(size_t b=0; b<bu.size(); b++){
      switch (bu[b].binvalue) {
        case binX:
        iscube = true;
        gridAxis1[0] = bu[b].min;
        gridAxis1[1] = bu[b].max;
        gridAxis1[2] = bu[b].bins();
        break;

        case binY:
        iscube = true;
        gridAxis2[0] = bu[b].min;
        gridAxis2[1] = bu[b].max;
        gridAxis2[2] = bu[b].bins();
        break;

        case binR:
        iscylinder = true;
        gridAxis1[0] = bu[b].min;
        gridAxis1[1] = bu[b].max;
        gridAxis1[2] = bu[b].bins();
        break;

        case binPhi:
        iscylinder = true;
        gridAxis2[0] = bu[b].min;
        gridAxis2[1] = bu[b].max;
        gridAxis2[2] = bu[b].bins();
        break;

        case binZ:
        gridAxis3[0] = bu[b].min;
        gridAxis3[1] = bu[b].max;
        gridAxis3[2] = bu[b].bins();
        break;

        case binRPhi:
        case binEta:
        case binH:
        case binMag:
        throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
        break;
      }
      if( !(iscylinder || iscube) || (iscylinder && iscube)){
        throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
      }
      if(iscylinder){
        transfoGlobalToLocal = [](Acts::Vector3D pos)->Acts::Vector3D{return{sqrt(pos.x()*pos.x()+pos.y()*pos.y()),
                                                      atan2(pos.y(),pos.x()), pos.z()}; };
        mapMaterial = mapMaterialCylinder;
      }
      if(iscube){
        transfoGlobalToLocal = [](Acts::Vector3D pos)->Acts::Vector3D{return{pos.x(), pos.y(), pos.z()};};
        mapMaterial = mapMaterial3D;
      }
    }
    MaterialGrid3D Grid = createMaterialGrid(gridAxis1, gridAxis2, gridAxis3, recMaterial.second, mapMaterial);
    // updateMapOverfow(Grid, mState.materialBin[recMaterial.first]);
    MaterialMapper<MaterialGrid3D> matMap(transfoGlobalToLocal, Grid);
    // for (size_t i = 0; i < Grid.numLocalBins()[0]+2; i++) {
    //   for (size_t j = 0; j < Grid.numLocalBins()[1]+2; j++) {
    //     for (size_t k = 0; k < Grid.numLocalBins()[2]+2; k++) {
    //       std::cout << "bin : " << i << " " << j << " " << k << std::endl;
    //       std::cout << "mat : " << Grid.atLocalBins({{i,j,k}}) << std::endl;
    //     }
    //   }
    // }
    mState.volumeMaterial[recMaterial.first] = std::make_unique<InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>>(std::move(matMap),std::move(mState.materialBin[recMaterial.first]));
  }
}


void Acts::VolumeMaterialMapper::mapMaterialTrack(
    State& mState, RecordedMaterialTrack& mTrack) const {
  // Neutral curvilinear parameters
  NeutralCurvilinearParameters start(std::nullopt, mTrack.first.first,
                                     mTrack.first.second, 0.);

  // Prepare Action list and abort list
  using DebugOutput = detail::DebugOutputActor;
  using MaterialVolumeCollector = VolumeCollector<MaterialVolume>;
  using ActionList = ActionList<MaterialVolumeCollector, DebugOutput>;
  using AbortList = AbortList<detail::EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options(mState.geoContext,
                                                   mState.magFieldContext);
  options.debug = m_cfg.mapperDebugOutput;

  // Now collect the material layers by using the straight line propagator
  const auto& result = m_propagator.propagate(start, options).value();
  auto mcResult = result.get<MaterialVolumeCollector::result_type>();
  // Massive screen output
  if (m_cfg.mapperDebugOutput) {
    auto debugOutput = result.get<DebugOutput::result_type>();
    ACTS_VERBOSE("Debug propagation output.");
    ACTS_VERBOSE(debugOutput.debugString);
  }

  auto mappingVolumes = mcResult.collected;

  // Retrieve the recorded material from the recorded material track
  auto& rMaterial = mTrack.second.materialInteractions;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material steps to map.")

  // These should be mapped onto the mapping surfaces found
  ACTS_VERBOSE("Found     " << mappingVolumes.size()
                            << " mapping volumes for this track.");
  ACTS_VERBOSE("Mapping volumes are :")
  for (auto& mVolumes : mappingVolumes) {
    ACTS_VERBOSE(" - Volume : " << mVolumes.volume->geoID()
                                 << " at position = (" <<
                                 mVolumes.position.x()
                                 << ", " << mVolumes.position.y() << ", "
                                 << mVolumes.position.z() << ")");

    mappingVolumes.push_back(mVolumes);
  }
  auto rmIter = rMaterial.begin();
  auto volIter = mappingVolumes.begin();
  bool encounterVolume = false;
  double mappingStep = 1.;

  int volumeStep = 1;
  Acts::Vector3D extraPosition = {0,0,0};
  Acts::Vector3D extraDirection = {0,0,0};


  while (rmIter != rMaterial.end() && volIter != mappingVolumes.end()) {
    if(volIter != mappingVolumes.end() && encounterVolume==true && !volIter->volume->inside(rmIter->position)){
      encounterVolume=false;
      ++volIter;
    }

    if(volIter != mappingVolumes.end() && volIter->volume->inside(rmIter->position)){
      volumeStep = floor(rmIter->materialProperties.thickness()/mappingStep);
      auto properties = rmIter->materialProperties;
      float remainder = rmIter->materialProperties.thickness() - mappingStep*volumeStep;
      properties.scaleThickness(mappingStep/properties.thickness());
      mState.recordedMaterial[volIter->volume->geoID()].push_back(std::pair(properties,rmIter->position));
      for(int step=1;step<=volumeStep;step++){
        extraDirection = rmIter->direction;
        extraDirection = extraDirection*(mappingStep/extraDirection.norm());
        extraPosition = rmIter->position + step*extraDirection;
        if(step==volumeStep){
          properties.scaleThickness(remainder/properties.thickness());
        }
        mState.recordedMaterial[volIter->volume->geoID()].push_back(std::pair(properties,extraPosition));
      }
      encounterVolume=true;
    }
    ++rmIter;
  }
}
