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
			const Config& config)
		: m_volume(volume), m_frustum(nullptr) {
			ACTS_VERBOSE("Constructing FrustumNavigationPolicy for volume " << m_volume.volumeName());
                        std::vector<BoundingBox*> prims;
                        m_boxes.push_back(std::make_unique<BoundingBox>(m_volume.boundingBox()));
                        prims.push_back(m_boxes.back().get());
                        for(auto & vol : m_volume.volumes()) {
				m_boxes.push_back(std::make_unique<BoundingBox>(vol.boundingBox()));
                                prims.push_back(m_boxes.back().get());
                        }
                        m_topBox=std::unique_ptr<BoundingBox>(make_octree(m_boxes,prims,config.depth));
		}

	void FrustumNavigationPolicy::initializeCandidates(const NavigationArguments& args,
			AppendOnlyNavigationStream& stream,
	      		const Logger& logger) const {
		ACTS_VERBOSE("FrustumNavigationPolicy Candidates initialization for volume"<< m_volume.volumeName());
                // Create the frustum if not present
		Frustum3 frustum(args.position,args.direction, M_PI / 4);
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
		const BoundingBox* topBoxCopy = m_topBox.get();
		while (topBoxCopy != nullptr) {
    			//if (topBoxCopy->intersect(*m_frustum.get())) {
			if(topBoxCopy->intersect(frustum)){
				if (topBoxCopy->hasEntity()) {
					const TrackingVolume* tvol=dynamic_cast<const TrackingVolume*>(topBoxCopy->entity());
					ACTS_VERBOSE("get portals from volume "<<tvol->volumeName());
					const auto& portals = tvol->portals();
					for(const auto& portal : portals){
						stream.addPortalCandidate(portal);
					}
					ACTS_VERBOSE("reset copy");
					topBoxCopy = topBoxCopy->getSkip();
				} else {
					ACTS_VERBOSE("no entity");
					topBoxCopy = topBoxCopy->getLeftChild();
				}
			} else {
				ACTS_VERBOSE("no intersection");
				topBoxCopy = topBoxCopy->getSkip();
			}
		}
		ACTS_VERBOSE("done with while loop");
	}

	void FrustumNavigationPolicy::connect(NavigationDelegate& delegate) const {
		connectDefault<FrustumNavigationPolicy>(delegate);
	}
}

