// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Frustum.hpp"

namespace Acts::Experimental {

	/// A navigation policy that uses frustum-octree intersections to find portals
	
	class FrustumNavigationPolicy : public INavigationPolicy {

		public:
			using BoundingBox = AxisAlignedBoundingBox<Volume, double, 3>;
			using Frustum3 = Frustum<double, 3, 3>;
			struct Config {
				/// The octree depth
				int depth=3;
			};
			/// Main constructor, which takes the top-level volume and builds the octree
			/// @param gctx The geometrycontext object
			/// @param volume The tracking volume used to build the octree
			/// @param config The configuration of the Navigation Policy
			/// @param logger A logging instance
			explicit FrustumNavigationPolicy(const GeometryContext& gctx, 
					const TrackingVolume& volume, 
					const Logger& logger,
					const Config& config)
			: m_volume(volume) {
				ACTS_VERBOSE("Constructing FrustumNavigationPolicy for volume " << m_volume.volumeName());
				std::vector<BoundingBox*> prims;
				m_boxes.push_back(std::make_unique<BoundingBox>(m_volume.boundingBox()));
				prims.push_back(m_boxes.back().get());
				for(auto & vol : m_volume.volumes()) {
					m_boxes.push_back(std::make_unique<BoundingBox>(vol.boundingBox()));
					prims.push_back(m_boxes.back().get());
				}
				m_topBox=make_octree(m_boxes,prims,config.depth);
			}

			/// Update the navigation state
  			/// @param args The navigation arguments
  			/// @param stream The navigation stream to update
  			/// @param logger The logger
  			void initializeCandidates(const NavigationArguments& args,
                            		AppendOnlyNavigationStream& stream,
                            		const Logger& logger) const {
    				ACTS_VERBOSE("FrustumNavigationPolicy Candidates initialization for volume"<< m_volume.volumeName());
				// Create the frustum if not present
				if(args.frustum == nullptr){
      					args.frustum=std::make_shared<Frustum3>(args.position,args.direction, M_PI / 4);
				}
    				// Check if we leave the frustum, reset candidates if so
    				const auto& normals = frustum->normals();
    				const bool outside = std::any_of(
						normals.begin(), normals.end(), [&args](const auto& normal) {
							return (args.position - frustum->origin()).dot(normal) >= 0;
						});
    				if (outside) {
					args.frustum=std::make_shared<Frustum3>(args.position,args.direction, M_PI / 4);
					//to be replaced by the new version of this when available
      					//args.surfaceCandidates.clear();
				}
				std::vector<const Portal*> portalCandidates = {};
				auto topBoxCopy = m_topBox;
    				while (topBoxCopy != nullptr) {
      					if (topBoxCopy->intersect(frustum)) {
        					if (topBoxCopy->hasEntity()) {
							const auto& portals = topBoxCopy->entity()->portals();
							for(const auto& portal : portals){
								stream.addPortalCandidate(portal);
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

			/// Connect this policy with a navigation delegate
			/// @param delegate The navigation delegate to connect to
  			void connect(NavigationDelegate& delegate) const override {
				connectDefault<FrustumNavigationPolicy>(delegate);
			}

		private:
			// The tracking volume
  			const TrackingVolume& m_volume;

			// The vector of unique pointers of bounding boxes
			std::vector<std::unique_ptr<BoundingBox> > m_boxes;

			// The top-level bounding box with the octree
			BoundingBox* m_topBox;
	};

	static_assert(NavigationPolicyConcept<FrustumNavigationPolicy>);
} //namespace Acts::Experimental
