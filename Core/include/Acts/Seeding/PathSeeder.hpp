// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/ISourceLinkGrid.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

/// @brief Seeding algorigthm that extracts 
/// the IP parameters and sorts the source links
/// into possible track candidates
///
/// @tparam axis_t Type of the axis to bin 
/// the source links
///
template <typename axis_t>
class PathSeeder {
    public:
        struct Seed {
            ActsScalar ipP;
            Vector3 ipDir;
            Vector3 ipVertex;
            std::vector<SourceLink> sourceLinks;
        };

        using SourceLinkCalibrator =
            Delegate<Vector3(
                const GeometryContext&,
                const SourceLink&)>;

        using PathWidthLookup = 
            Delegate<std::pair<ActsScalar, ActsScalar>(
                const GeometryContext&, 
                const GeometryIdentifier&)>;

        using IntersectionLookup = 
            Delegate<std::vector<std::pair<GeometryIdentifier,Vector3>>(
                const GeometryContext&, 
                const Vector3&, 
                const Vector3&,
                const ActsScalar&)>;

        using TrackEstimator = 
            Delegate<std::tuple<ActsScalar, Vector3, Vector3, Vector3>(
                const GeometryContext&, 
                const Vector3&)>;

        /// @brief The nested configuration struct
        struct Config {
            /// Binned SourceLink provider
            std::shared_ptr<ISourceLinkGrid<axis_t>> sourceLinkGrid;
            /// Parameters estimator
            TrackEstimator trackEstimator;
            /// Surface accessor
            SourceLinkCalibrator sourceLinkCalibrator;
            /// Intersection finder
            IntersectionLookup intersectionFinder;
            /// Path width provider
            PathWidthLookup pathWidthProvider;
            /// First layer extent
            Extent firstLayerExtent;
            /// Direction of the telescope extent
            BinningValue orientation = binX;
        };

        /// @brief Constructor
        PathSeeder(const Config& config)
            : m_cfg(std::move(config)) {};

        /// @brief Destructor
        ~PathSeeder() = default;


        /// @brief Extract the IP parameters and 
        /// sort the source links into the seeds
        ///
        /// @param gctx The geometry context
        /// @param sourceLinks The source links to seed
        ///
        /// @return The vector of seeds
        std::vector<Seed>
        getSeeds(
            const GeometryContext& gctx, 
            const std::vector<SourceLink>& sourceLinks) const {
                // Sort the source links into the grid
                m_cfg.sourceLinkGrid->initialize(gctx, sourceLinks);

                // Get plane of the telescope
                // sensitive surfaces
                BinningValue bin0, bin1;
                if (m_cfg.orientation == binX) {
                    bin0 = binY;
                    bin1 = binZ;
                } else if (m_cfg.orientation == binY) {
                    bin0 = binX;
                    bin1 = binZ;
                } else {
                    bin0 = binX;
                    bin1 = binY;
                }

                // Create the seeds
                std::vector<Seed> seeds;
                for (const auto& sl : sourceLinks) {
                    Vector3 globalPos = m_cfg.sourceLinkCalibrator(
                        gctx, sl);

                    // Check if the hit is in the 
                    // first tracking layer
                    if (!m_cfg.firstLayerExtent.contains(globalPos)) {
                        continue;
                    }

                    // Get the IP parameters
                    auto [ipP, ipVertex, ipDir, flDir] = 
                        m_cfg.trackEstimator(
                            gctx, globalPos);

                    // Intersect with the surfaces
                    std::vector<std::pair<GeometryIdentifier,Vector3>>
                        intersections =
                            m_cfg.intersectionFinder(
                                gctx, globalPos, flDir, ipP);
            
                    // Continue if no intersections
                    if (intersections.empty()) {
                        continue;
                    }
                    // Vector to store the source links
                    std::vector<SourceLink> seedSourceLinks;

                    // Store the IP parameters
                    Seed seed;
                    seed.ipP = ipP;
                    seed.ipDir = ipDir;
                    seed.ipVertex = ipVertex;

                    // Store the pivot source link
                    seedSourceLinks.push_back(sl);
                
                    // Iterate over the intersections 
                    // and get the source links
                    // in the subsequent layers
                    for (auto& [geoId,refPoint] : intersections) {
                        // Get the path width
                        auto [pathWidth0, pathWidth1] = 
                            m_cfg.pathWidthProvider(
                                gctx, 
                                geoId);

                        // Get the bounds of the path
                        ActsScalar top0 = 
                            refPoint[bin0] + pathWidth0;
                        ActsScalar bot0 = 
                            refPoint[bin0] - pathWidth0;
                        ActsScalar top1 = 
                            refPoint[bin1] + pathWidth1;
                        ActsScalar bot1 = 
                            refPoint[bin1] - pathWidth1;
                
                        // Get the lookup table for the source links
                        auto grid = m_cfg.sourceLinkGrid->getSourceLinkTable(
                            geoId);
                
                        // Get the range of bins to search for source links
                        auto botLeftBin = grid.localBinsFromPosition(
                            Vector2(bot0, bot1));
                        auto topRightBin = grid.localBinsFromPosition(
                            Vector2(top0, top1));
                
                        // Get the source links from the lookup table
                        // by iterating over the bin ranges
                        auto currentBin = botLeftBin;
                        while (currentBin.at(1) <= topRightBin.at(1)) {
                            while (currentBin.at(0) <= topRightBin.at(0)) {
                                auto sourceLinksToAdd = 
                                    grid.atLocalBins(currentBin);

                                seedSourceLinks.insert(
                                    seedSourceLinks.end(),
                                    sourceLinksToAdd.begin(),
                                    sourceLinksToAdd.end());
                                currentBin.at(0)++;
                            }
                            currentBin.at(1)++;
                            currentBin.at(0) = botLeftBin.at(0);
                        }
                    }

                    // Add the source links to the seed
                    seed.sourceLinks = seedSourceLinks;

                    // Add the seed to the list
                    seeds.push_back(seed);
                }
                return seeds;
        };


    private:
        Config m_cfg;
};

} // namespace Acts
