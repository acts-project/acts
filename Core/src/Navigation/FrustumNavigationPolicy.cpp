// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/FrustumNavigationPolicy.hpp"

namespace Acts::Experimental {

FrustumNavigationPolicy::FrustumNavigationPolicy(const GeometryContext& gctx,
                                                 const TrackingVolume& volume,
                                                 const Logger& logger,
                                                 const Config& config) {
  ACTS_DEBUG("Constructing FrustumNavigationPolicy for volume "
             << volume.volumeName());
  m_name = volume.volumeName();
  std::vector<BoundingBox*> prims;
  m_boxes.push_back(std::make_unique<BoundingBox>(volume.boundingBox(gctx)));
  prims.push_back(m_boxes.back().get());
  for (auto& vol : volume.volumes()) {
    ACTS_DEBUG("add volume " << vol.volumeName()
                             << " to list of bounding boxes");
    m_boxes.push_back(std::make_unique<BoundingBox>(vol.boundingBox(gctx)));
    prims.push_back(m_boxes.back().get());
  }
  m_topBox = make_octree(m_boxes, prims, config.depth);
}

void FrustumNavigationPolicy::initializeCandidates(
    const GeometryContext& gctx, const NavigationArguments&,
    NavigationPolicyState& state, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
  ACTS_DEBUG("FrustumNavigationPolicy Candidates initialization for volume "
             << m_name);
  auto& s = state.as<State>();
  ACTS_DEBUG("Frustum origin " << s.frustum.origin() << ", frustum dir "
                               << s.frustum.dir());
  const BoundingBox* topBoxCopy = m_topBox;
  while (topBoxCopy != nullptr) {
    if (topBoxCopy->intersect(s.frustum)) {
      if (topBoxCopy->hasEntity()) {
        const TrackingVolume* tvol =
            dynamic_cast<const TrackingVolume*>(topBoxCopy->entity());
        ACTS_DEBUG("get portals from volume " << tvol->volumeName());
        const auto& portals = tvol->portals();
        for (const auto& portal : portals) {
	  //To avoid including unnecessary portals, skip any portals from the top-level volume to its child volumes
	  //If we want these, they will be added when the intersection reaches the child volume
	  //Only add portals from the top-level volume to volumes it doesn't contain
          if (tvol->volumeName().compare(m_name) ==
              0) {
            Acts::Result<const TrackingVolume*> pvolr =
                portal.resolveVolume(gctx, s.frustum.origin(), s.frustum.dir());
            if (pvolr.ok()) {
              const TrackingVolume* pvol = *pvolr;
              if (tvol->inside(gctx, pvol->center(gctx))) {
                continue;
              } else {
                stream.addPortalCandidate(portal);
              }
            } else {
              ACTS_DEBUG("unable to resolve portal volume"); }
          } else {
            stream.addPortalCandidate(portal);}
        }
        topBoxCopy = topBoxCopy->getSkip();
      } else {
        topBoxCopy = topBoxCopy->getLeftChild();
      }
    } else {
      topBoxCopy = topBoxCopy->getSkip();
    }
  }
}

void FrustumNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<FrustumNavigationPolicy>(delegate);
}

bool FrustumNavigationPolicy::isValid(const GeometryContext&,
                                      const NavigationArguments& args,
                                      NavigationPolicyState& state,
                                      const Logger& logger) const {
  // Check if we leave the frustum, reset candidates if so
  auto& s = state.as<State>();
  const auto& difference = args.position - s.frustum.origin();
  const auto& normals = s.frustum.normals();
  auto it_start = std::next(normals, 1);
  const bool outside = std::any_of(
      it_start, normals.end(),
      [&difference](const auto& normal) { return difference.dot(normal) < 0; });
  if (outside) {
    ACTS_DEBUG("FrustumNavigationPolicy: outside frustum");
    return false;
  } else {
    ACTS_DEBUG("FrustumNavigationPolicy: frustum still ok");
    return true;
  }
}

void FrustumNavigationPolicy::createState(
    const GeometryContext&, const NavigationArguments& args,
    NavigationPolicyStateManager& stateManager, const Logger& logger) const {
  ACTS_DEBUG("create FrustumNavigationPolicy state");
  auto& s = stateManager.pushState<State>();
  s.frustum = Frustum3(args.position, args.direction, std::numbers::pi / 2);
}

void FrustumNavigationPolicy::popState(
    NavigationPolicyStateManager& stateManager, const Logger& logger) const {
  ACTS_DEBUG("remove FrustumNavigationPolicy state");
  stateManager.popState();
}
}  // namespace Acts::Experimental
