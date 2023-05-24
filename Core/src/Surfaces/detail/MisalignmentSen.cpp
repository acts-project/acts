#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <iostream>
#include <unordered_map>
#include <array>

// using array to define misalignment parameters (fixed size)
using SensorMisalignment = std::array<double, 6>;

namespace ActsAlignment {
namespace detail {

using namespace Acts;

// ...


struct TrackAlignmentState {
  
  // The misalignment parameters for each aligned surface
  std::unordered_map<const Surface*, SensorMisalignment> misalignments;
};

// Define a function to obtain misalignment parameters for a surface
SensorMisalignment getMisalignmentParametersForSurface(const Surface* surface) {
  // Implementation of the function to retrieve misalignment parameters for a surface
  // ...
  return SensorMisalignment{}; 

// @TO DO: how to add changes regarding the misalignment
template <typename traj_t, typename parameters_t = BoundTrackParameters>
TrackAlignmentState trackAlignmentState(
    const GeometryContext& gctx, const Acts::MultiTrajectory<traj_t>& multiTraj,
    const size_t& entryIndex,
    const std::pair<ActsDynamicMatrix, std::unordered_map<size_t, size_t>>&
        globalTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces,
    const AlignmentMask& alignMask) {
  // ...

  TrackAlignmentState alignState; // Define and initialize alignState variable

  // Loop over the measurement states to fill those alignment matrices
  // This is done in reverse order - important!
  size_t iMeasurement = alignState.measurementDim;
  size_t iParams = alignState.trackParametersDim;
  size_t iSurface = nAlignSurfaces;

  for (const auto& [rowStateIndex, isAlignable] : measurementStates) {

    // (d) Get the derivative of alignment parameters w.r.t. measurement
    // or residual
    if (isAlignable) {
      iSurface -= 1;
      const auto surface = &state.referenceSurface();
      alignState.alignedSurfaces.at(surface).second = iSurface;

      // For each surface, we're trying to store misalignment parameters associated to this surface
      SensorMisalignment misalignment = getMisalignmentParametersForSurface(surface);

      // Apply the misalignment to the surface
      surface->applyMisalignment(misalignment[0], misalignment[1],
                                 misalignment[2], misalignment[3],
                                 misalignment[4], misalignment[5]);

      // alignment state = 'place' where we can store misalignment parameters
      alignState.misalignments[surface] = misalignment;
    }
  }

  return alignState;
}

}  
}  

using namespace ActsAlignment::detail;


// for specific surface, return the misalign parameters 
SensorMisalignment getMisalignmentParametersForSurface(const Surface* surface) {
  double DeltaX = getDeltaXForSurface(surface);
  double DeltaY = getDeltaYForSurface(surface);
  double DeltaZ = getDeltaZForSurface(surface);
  double DeltaPhi = getDeltaPhiForSurface(surface);
  double DeltaTheta = getDeltaThetaForSurface(surface);
  double DeltaPsi = getDeltaPsiForSurface(surface);

  // Then we can create the SensorMisalignment object that will be initialise with this set of paramteres 
  SensorMisalignment misalignment = {DeltaX, DeltaY, DeltaZ,
                                      DeltaPhi, DeltaTheta, DeltaPsi};
  return misalignment;
}

int main() {

  //  Setting up - initialisation and definition of variables that we will need
  GeometryContext geometryContext;
  MultiTrajectory<traj_t> multiTrajectory;
  size_t entryIndex = 0;
  std::pair<ActsDynamicMatrix, std::unordered_map<size_t, size_t>> globalTrackParamsCov;
  std::unordered_map<const Surface*, size_t> indexedAlignSurfaces;
  AlignmentMask alignMask;

  // Call the trackAlignmentState function with your input data
  TrackAlignmentState alignmentState = trackAlignmentState(
      geometryContext, multiTrajectory, entryIndex, globalTrackParamsCov,
      indexedAlignSurfaces, alignMask);

  for (const auto& pair : indexedAlignSurfaces) {
    const Surface* surface = pair.first;
    SensorMisalignment misalignment = alignmentState.misalignments[surface];
  
    std::cout << "Misalignment parameters for the surface:" << std::endl;
    std::cout << "Index of a surface: " << pair.second << std::endl;
    std::cout << "Delta X: " << misalignment[0] << std::endl;
    std::cout << "Delta Y: " << misalignment[1] << std::endl;
    std::cout << "Delta Z: " << misalignment[2] << std::endl;
    std::cout << "Delta Phi: " << misalignment[3] << std::endl;
    std::cout << "Delta Theta: " << misalignment[4] << std::endl;
    std::cout << "Delta Psi: " << misalignment[5] << std::endl;
  }

  // This should work only if we're assuming that the index of the layer is known
  size_t layerIndex = 0;  

  // in indexedAlignSurfaces, try to find a surface that belongs to the specific layer
  const Surface* layerSurface = nullptr;
  for (const auto& pair : indexedAlignSurfaces) {
    const Surface* surface = pair.first;
    // Assuming the layer index is stored in the surface object
    if (surface->layerIndex == layerIndex) {
      layerSurface = surface;
      break;
    }
  }

  // If this layer is found, we can get the misalignment parameters (layer surface)
  if (layerSurface != nullptr) {
    SensorMisalignment misalignment = alignmentState.misalignments[layerSurface];
    std::cout << "Layer surface: - list of misalignment parameters:" << std::endl;
    std::cout << "Delta X: " << misalignment[0] << std::endl;
    std::cout << "Delta Y: " << misalignment[1] << std::endl;
    std::cout << "Delta Z: " << misalignment[2] << std::endl;
    std::cout << "Delta Phi: " << misalignment[3] << std::endl;
    std::cout << "Delta Theta: " << misalignment[4] << std::endl;
    std::cout << "Delta Psi: " << misalignment[5] << std::endl;
  } else {
    std::cout << "Layer surface can't be found." << std::endl;
  }

  return 0;
}
