// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/FrustumNavigationPolicy.hpp"

#include "Acts/Geometry/TrivialPortalLink.hpp"

namespace Acts::Experimental {

namespace {

/// Helper function  to check whether the portal sits on the original boundary 
/// surface of the volume
/// @param link: Pointer to the portallink to check
/// @param volume: Reference to the volume which boundaries need to be found 
bool leadsOutside(const PortalLinkBase *link, const TrackingVolume &volume) {
  const auto *trivialLink = dynamic_cast<const TrivialPortalLink *>(link);
  if (trivialLink == nullptr) {
    return false;
  }

  const TrackingVolume *volumeLink = &trivialLink->volume();
  if (volumeLink == &volume || volumeLink->motherVolume() == &volume) {
    return false;
  }
  return true;
}

}  // namespace

FrustumNavigationPolicy::FrustumNavigationPolicy(const GeometryContext &gctx,
                                                 const TrackingVolume &volume,
                                                 const Logger &logger,
                                                 const Config &config) {
  ACTS_VERBOSE("Constructing FrustumNavigationPolicy for volume "
               << volume.volumeName());
  m_id = volume.geometryId();
  std::vector<BoundingBox *> prims;
  m_boxes.push_back(std::make_unique<BoundingBox>(volume.boundingBox(gctx)));
  prims.push_back(m_boxes.back().get());
  for (auto &vol : volume.volumes()) {
    ACTS_VERBOSE("add volume " << vol.volumeName()
                               << " to list of bounding boxes");
    m_boxes.push_back(std::make_unique<BoundingBox>(vol.boundingBox(gctx)));
    prims.push_back(m_boxes.back().get());
  }
  m_topBox =
      Acts::BoundingBoxHierarchy::makeOctree(m_boxes, prims, config.depth);

  for (const auto &portal : volume.portals()) {
    if (leadsOutside(portal.getLink(Direction::AlongNormal()), volume) ||
        leadsOutside(portal.getLink(Direction::OppositeNormal()), volume)) {
      m_cachedPortals.push_back(&portal);
    }
  }
}

void FrustumNavigationPolicy::initializeCandidates(
    const GeometryContext & /*gctx*/, const NavigationArguments &args,
    NavigationPolicyState &state, AppendOnlyNavigationStream &stream,
    const Logger &logger) const {
  ACTS_VERBOSE("FrustumNavigationPolicy Candidates initialization for volume "
               << m_id);
  auto &s = state.as<State>();
  s.frustum = Frustum3(args.position, args.direction, s.openingAngle);

  ACTS_VERBOSE("Frustum origin " << s.frustum.origin() << ", frustum dir "
                                 << s.frustum.dir());
  Acts::BoundingBoxHierarchy::visitIntersecting(
      s.frustum, m_topBox, [this, &stream, &logger](const Volume &entity) {
        const TrackingVolume *tvol =
            dynamic_cast<const TrackingVolume *>(&entity);
        ACTS_VERBOSE("Get portals from volume " << tvol->volumeName());
        if (tvol->geometryId() == m_id) {
          for (const Portal *portal : m_cachedPortals) {
            stream.addPortalCandidate(*portal);
          }
        } else {
          for (const auto &portal : tvol->portals()) {
            stream.addPortalCandidate(portal);
          }
        }
      });
  ACTS_VERBOSE(
      "FrustumNavigationPolicy Candidates initialization done for volume "
      << this->m_id);
}

void FrustumNavigationPolicy::connect(NavigationDelegate &delegate) const {
  connectDefault<FrustumNavigationPolicy>(delegate);
}

bool FrustumNavigationPolicy::isValid(const GeometryContext & /*gctx*/,
                                      const NavigationArguments &args,
                                      NavigationPolicyState &state,
                                      const Logger &logger) const {
  // Check if we leave the frustum, reset candidates if so
  auto &s = state.as<State>();
  double costheta =
      args.position.normalized().dot(s.frustum.origin().normalized());
  // Compare to half the frustum opening angle
  if (costheta < std::cos(s.openingAngle / 2)) {
    ACTS_DEBUG("FrustumNavigationPolicy: outside frustum");
    return false;
  } else {
    ACTS_DEBUG("FrustumNavigationPolicy: frustum still ok");
    return true;
  }
}

void FrustumNavigationPolicy::createState(
    const GeometryContext & /*gctx*/, const NavigationArguments &args,
    NavigationPolicyStateManager &stateManager, const Logger &logger) const {
  ACTS_DEBUG("create FrustumNavigationPolicy state");
  auto &s = stateManager.pushState<State>();
  s.frustum = Frustum3(args.position, args.direction, std::numbers::pi / 4);
  s.openingAngle = std::numbers::pi / 4;
}

void FrustumNavigationPolicy::popState(
    NavigationPolicyStateManager &stateManager, const Logger &logger) const {
  ACTS_DEBUG("remove FrustumNavigationPolicy state");
  stateManager.popState();
}
}  // namespace Acts::Experimental
