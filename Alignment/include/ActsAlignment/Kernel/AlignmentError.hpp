// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// acts/Alignment/include/ActsAlignment/Kernel/AlignmentError.hpp

#pragma once

#include <iostream>
#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace ActsAlignment {
// This is the custom error code enum
enum class AlignmentError {
  NoAlignmentDofOnTrack = 1,
  AlignmentParametersUpdateFailure = 2,
  ConvergeFailure = 3
};

namespace detail {
// Define a custom error code category derived from std::error_category
class AlignmentErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "AlignmentError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<AlignmentError>(c)) {
      case AlignmentError::NoAlignmentDofOnTrack:
        return "No alignment parameters on the track";
      case AlignmentError::AlignmentParametersUpdateFailure:
        return "Update to alignment parameters failure";
      case AlignmentError::ConvergeFailure:
        return "The alignment is not converged";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::AlignmentErrorCategory& AlignmentErrorCategory() {
  static detail::AlignmentErrorCategory c;
  return c;
}

inline std::error_code make_error_code(ActsAlignment::AlignmentError e) {
  return {static_cast<int>(e), ActsAlignment::AlignmentErrorCategory()};
}

struct MisalignmentParameters {
  Acts::Translation3 translation;
  Acts::RotationMatrix3 rotation;
};

template <typename fitter_t>
template <typename trajectory_container_t,
          typename start_parameters_container_t, typename fit_options_t>
Acts::Result<ActsAlignment::AlignmentResult>
ActsAlignment::Alignment<fitter_t>::align(
    const trajectory_container_t& trajectoryCollection,
    const start_parameters_container_t& startParametersCollection,
    const ActsAlignment::AlignmentOptions<fit_options_t>& alignOptions) const {
  // Construct an AlignmentResult object
  AlignmentResult alignResult;

  // Assign index to the alignable surface
  for (unsigned int iDetElement = 0;
       iDetElement < alignOptions.alignedDetElements.size(); iDetElement++) {
    alignResult.idxedAlignSurfaces.emplace(
        &alignOptions.alignedDetElements.at(iDetElement)->surface(),
        iDetElement);
  }
  ACTS_VERBOSE("There are " << alignResult.idxedAlignSurfaces.size()
                            << " detector elements to be aligned");

  // Start the iteration to minimize the chi2
  bool converged = false;
  bool alignmentParametersUpdated = false;
  std::queue<double> recentChi2ONdf;
  ACTS_INFO("Max number of iterations: " << alignOptions.maxIterations);
  for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
    // Perform the fit to the trajectories and update alignment parameters
    // Initialize the alignment mask (all dof in default)
    AlignmentMask alignMask = AlignmentMask::All;
    // Set the alignment mask
    auto iter_it = alignOptions.iterationState.find(iIter);
    if (iter_it != alignOptions.iterationState.end()) {
      alignMask = iter_it->second;
    }
    // Calculate the alignment parameters delta etc.
    calculateAlignmentParameters(
        trajectoryCollection, startParametersCollection,
        alignOptions.fitOptions, alignResult, alignMask);
    // Screen out the information
    ACTS_INFO("iIter = " << iIter << ", total chi2 = " << alignResult.chi2
                         << ", total measurementDim = "
                         << alignResult.measurementDim
                         << " and average chi2/ndf = "
                         << alignResult.averageChi2ONdf);
    // Check if it has converged against the provided precision
    // 1. either the delta average chi2/ndf in the last few
    // iterations is within tolerance
    if (recentChi2ONdf.size() >=
        alignOptions.deltaAverageChi2ONdfCutOff.first) {
      if (std::abs(recentChi2ONdf.front() - alignResult.averageChi2ONdf) <=
          alignOptions.deltaAverageChi2ONdfCutOff.second) {
        ACTS_INFO(
            "Alignment has converged with change of chi2/ndf < "
            << alignOptions.deltaAverageChi2ONdfCutOff.second << " in the last "
            << alignOptions.deltaAverageChi2ONdfCutOff.first << " iterations"
            << " after " << iIter << " iteration(s)");
        converged = true;
        break;
      }
      recentChi2ONdf.pop();
    }
    // 2. or the average chi2/ndf (is this correct?)
    if (alignResult.averageChi2ONdf <= alignOptions.averageChi2ONdfCutOff) {
      ACTS_INFO("Alignment has converged with average chi2/ndf < "
                << alignOptions.averageChi2ONdfCutOff << " after " << iIter
                << " iteration(s)");
      converged = true;
      break;
    }
    // Remove the first element
    // Store the result in the queue
    recentChi2ONdf.push(alignResult.averageChi2ONdf);

    ACTS_INFO("The solved delta of alignmentParameters = \n "
              << alignResult.deltaAlignmentParameters);
    // Not converged yet, update the detector element alignment parameters
    auto updateRes = updateAlignmentParameters(
        alignOptions.fitOptions.geoContext, alignOptions.alignedDetElements,
        alignOptions.alignedTransformUpdater, alignResult);
    if (!updateRes.ok()) {
      ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
      return updateRes.error();
    }
    alignmentParametersUpdated = true;
  }  // end of all iterations

  // Alignment failure if not converged
  if (!converged) {
    ACTS_ERROR("Alignment is not converged.");
    alignResult.result = AlignmentError::ConvergeFailure;
  }

  // Apply correlated misalignments to the sensors
  if (alignmentParametersUpdated) {
    for (const auto& [surface, index] : alignResult.idxedAlignSurfaces) {
      const MisalignmentParameters& misalignmentParams =
          alignOptions.misalignmentParameters.at(index);

      // Calculate the correlated misalignment parameters
      Acts::AlignmentVector correlatedMisalignmentParams =
          Acts::AlignmentVector::Zero();
      correlatedMisalignmentParams.segment<3>(Acts::eAlignmentCenter0) =
          misalignmentParams.translation;
      correlatedMisalignmentParams.segment<3>(Acts::eAlignmentRotation0) =
          misalignmentParams.rotation.eulerAngles(2, 1, 0);


// Set the updated transform for the sensor
sensor->setTransform(updatedTransform);

ACTS_VERBOSE("Sensor with surface " << surface->geometryId()
                                    << " has aligned geometry position as below:");
ACTS_VERBOSE("Center (cenX, cenY, cenZ) = " << correlatedTranslation.transpose());
ACTS_VERBOSE("Euler angles (rotZ, rotY, rotX) = " << correlatedMisalignmentParams.transpose());
ACTS_VERBOSE("Rotation matrix = \n" << correlatedRotation);
} // end of for loop for sensors
} else {
ACTS_DEBUG("Alignment parameters are not updated.");
}

return alignResult;
} // end of align function
} // namespace ActsAlignment
