// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

namespace Acts {



void MultiLayerNavigationPolicy::initializeCandidates(
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
        // Custom implementation
    ACTS_VERBOSE("MultiLayerNavPol (volume=" << volume().volumeName() << ")");
  
    //create an artificial path inside the chamber to intersect the grid surfaces
    //express the path in the local frame
    const Transform3& transform = volume().transform();
    SurfaceArray* surfArray = surfaceArray();

    const Acts::Vector3 locPosition = transform*
        args.position;
    const Acts::Vector3 locDirection = transform.linear()*
        args.direction;
    
    auto step =  std::sqrt(std::pow(surfArray->binWidth(1u)[0], 2) +
    std::pow(surfArray->binWidth(1u)[1], 2));
    auto nSteps = surfArray->getAxes()[1]->getNBins();

    //generate the projected path on the grid - positions in the local frame
    std::vector<Vector3> path = generatePath(locPosition, locDirection, step, nSteps);

    //extract the sensitive sufrcaes from the grid that are along the path

    std::vector<const Surface*> surfCandidates = {};

    for(const auto& pos : path){
        //get the sensitive surfaces with the neighbors as configured from the bin extension
        const std::vector<const Surface*>& sensitiveSurfaces =
        surfArray->neighbors(pos);
        ACTS_VERBOSE("~> MultiLayer Navigation Policy reports" << sensitiveSurfaces.size()
                                             << " sensitive surfaces");
        surfCandidates.insert(surfCandidates.end(), sensitiveSurfaces.begin(), sensitiveSurfaces.end());
  
    }

    //fill the navigation stream with the surface candidates
    for(const auto *surf : surfCandidates){
        stream.addSurfaceCandidate(*surf, args.tolerance);

    }

}

void MultiLayerNavigationPolicy::connect(NavigationDelegate& delegate) const {
        // Custom implementation
        connectDefault<MultiLayerNavigationPolicy>(delegate);
}

std::vector<Vector3> MultiLayerNavigationPolicy::generatePath(
    const Vector3& startPosition, const Vector3& direction, double stepSize,
    std::size_t numberOfSteps) const {
        // Custom implementation
        std::vector<Vector3> path;
        path.reserve(numberOfSteps);
        for (std::size_t i = 0; i < numberOfSteps; ++i) {
            path.push_back(startPosition + i * stepSize * direction);
        }
        return path;

}
} // namespace Acts