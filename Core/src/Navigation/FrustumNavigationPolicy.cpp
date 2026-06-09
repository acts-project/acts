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
			const Config& config){
			ACTS_INFO("Constructing FrustumNavigationPolicy for volume " << volume.volumeName());
			m_name=volume.volumeName();
                        std::vector<BoundingBox*> prims;
                        m_boxes.push_back(std::make_unique<BoundingBox>(volume.boundingBox(volume.center(gctx))));
                        prims.push_back(m_boxes.back().get());
                        for(auto & vol : volume.volumes()) {
				ACTS_DEBUG("add volume "<<vol.volumeName()<<" to list of bounding boxes");
				BoundingBox bb=vol.boundingBox();
				m_boxes.push_back(std::make_unique<BoundingBox>(&vol,bb.min()+vol.center(gctx),bb.max()+vol.center(gctx)));
                                prims.push_back(m_boxes.back().get());
                        }
                        m_topBox=make_octree(m_boxes,prims,config.depth);
		}

	void FrustumNavigationPolicy::initializeCandidates(const GeometryContext& gctx,
			const NavigationArguments& args,
			NavigationPolicyState& state,
			AppendOnlyNavigationStream& stream,
	      		const Logger& logger) const {
		ACTS_INFO("FrustumNavigationPolicy Candidates initialization for volume "<< m_name);
		auto& s = state.as<State>();
		ACTS_DEBUG("Frustum origin "<<s.frustum.origin()<<", frustum dir "<<s.frustum.dir());
                // Recreate the frustum - needed?
		//s.frustum=Frustum3(args.position,args.direction, std::numbers::pi / 4);
		//ACTS_INFO("Frustum origin "<<s.frustum.origin()<<", frustum dir "<<s.frustum.dir());
                /*
		if(m_frustum.get() == nullptr){
			m_frustum.reset(new Frustum3(args.position,args.direction, M_PI / 4));
		}
		// Check if we leave the frustum, reset candidates if so
		const auto& difference=args.position-m_frustum->origin();
		const auto& normals = m_frustum->normals();
		const auto& difference=args.position-frustum->origin();
		const bool outside = std::any_of(
				normals.begin(), normals.end(), [&difference](const auto& normal) {
				return difference.dot(normal) >= 0;
				});
		if (outside) {
			m_frustum.reset(new Frustum3(args.position,args.direction, M_PI / 4));
			//to be replaced by the new version of this when available
			//args.surfaceCandidates.clear();
		}
		*/
		const BoundingBox* topBoxCopy = m_topBox;
		ACTS_DEBUG("search the octree");
		while (topBoxCopy != nullptr) {
			if(topBoxCopy->intersect(s.frustum)){
				if (topBoxCopy->hasEntity()) {
					const TrackingVolume* tvol=dynamic_cast<const TrackingVolume*>(topBoxCopy->entity());
					ACTS_DEBUG("get portals from volume "<<tvol->volumeName());
					const auto& portals = tvol->portals();
					ACTS_DEBUG("add "<<portals.size()<<" portal candidates");
					for(const auto& portal : portals){
						if(tvol->volumeName().compare(m_name)==0){ //only check for top-level volume
							Acts::Result<const TrackingVolume*> pvolr=portal.resolveVolume(gctx,s.frustum.origin(),s.frustum.dir());
							if(pvolr.ok()){
								const TrackingVolume* pvol=*pvolr;
								if(tvol->inside(gctx,pvol->center(gctx))){
									ACTS_DEBUG("skip portal to volume "<<pvol->volumeName());
									continue;
								}
								else{
									ACTS_DEBUG("add portal to volume "<<pvol->volumeName());
									stream.addPortalCandidate(portal);
								}
							}
							else ACTS_DEBUG("unable to resolve portal volume");
						}
						else stream.addPortalCandidate(portal);
						//stream.addPortalCandidate(portal);
					}
					ACTS_DEBUG("reset copy");
					topBoxCopy = topBoxCopy->getSkip();
				} else {
					ACTS_DEBUG("no entity");
					topBoxCopy = topBoxCopy->getLeftChild();
				}
			} else {
				if (topBoxCopy->hasEntity()) {
					const TrackingVolume* tvol=dynamic_cast<const TrackingVolume*>(topBoxCopy->entity());
					ACTS_DEBUG("no intersection with "<<tvol->volumeName());
				}
				topBoxCopy = topBoxCopy->getSkip();
			}
		}
		ACTS_DEBUG("done with while loop");
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
		const auto& difference=args.position-s.frustum.origin();
		//ACTS_DEBUG("FrustumNavigationPolicy: position is "<<args.position);
		//ACTS_DEBUG("FrustumNavigationPolicy: direction is "<<args.direction);
		//ACTS_DEBUG("FrustumNavigationPolicy: frustum origin is "<<s.frustum.origin());
		//for(auto& norm : s.frustum.normals()){
		//	ACTS_DEBUG("FrustumNavigationPolicy: frustum normal "<<norm);
		//	ACTS_DEBUG("FrustumNavigationPolicy: dot product="<<difference.dot(norm));
		//}
		const auto& normals = s.frustum.normals();
		auto it_start=normals.begin(); ++it_start;
		const bool outside = std::any_of(
				it_start, normals.end(), [&difference](const auto& normal) {
				return difference.dot(normal) < 0;
				});
		if (outside){
			ACTS_DEBUG("FrustumNavigationPolicy: outside frustum");
		       	return false;
		}
		else{
			ACTS_DEBUG("FrustumNavigationPolicy: frustum still ok");
		       	return true;
		}
	}

	void FrustumNavigationPolicy::createState(const GeometryContext&,
                         const NavigationArguments& args,
                         NavigationPolicyStateManager& stateManager,
                         const Logger& logger) const {
		ACTS_DEBUG("create FrustumNavigationPolicy state");
		auto& s = stateManager.pushState<State>();
		s.frustum=Frustum3(args.position,args.direction,std::numbers::pi / 2);
	}

	void FrustumNavigationPolicy::popState(NavigationPolicyStateManager& stateManager,
                      const Logger& logger) const {
		ACTS_DEBUG("remove FrustumNavigationPolicy state");
		stateManager.popState();
	}
}

