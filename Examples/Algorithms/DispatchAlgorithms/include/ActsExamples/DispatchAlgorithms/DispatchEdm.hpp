// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <vector>
#include "ActsExamples/EventData/Index.hpp"

namespace ActsExamples {

    struct DispatchMeasurements {
        // The index into the cluster collection
        std::vector<int> clusterIndices;

        // The geoId of the surface
        std::vector<int> clusterGeoIds;
        
        // Global information
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        // Local information
        std::vector<double> lx;
        std::vector<double> ly;
        std::vector<double> covLxx;
        std::vector<double> covLyy;

        /// @brief Reserve space for n measurements in all vectors
        /// @param n the size to reserve
        void reserve(std::size_t n) {
            clusterIndices.reserve(n);
            clusterGeoIds.reserve(n);
            x.reserve(n);
            y.reserve(n);
            z.reserve(n);
            lx.reserve(n);
            ly.reserve(n);
            covLxx.reserve(n);
            covLyy.reserve(n);
        }
        
    };

    struct DispatchParticles {
        std::vector<int> particlePdgs;
        std::vector<double> px;
        std::vector<double> py;
        std::vector<double> pz;
        std::vector<double> vx;
        std::vector<double> vy;
        std::vector<double> vz;

        /// @brief Reserve space for n particles in all vectors
        /// @param n the size to reserve
        void reserve(std::size_t n) {
            particlePdgs.reserve(n);
            px.reserve(n);
            py.reserve(n);
            pz.reserve(n);
            vx.reserve(n);
            vy.reserve(n);
            vz.reserve(n);
        }
    };

    using DispatchParticleMeasurementsMap = InverseMultimap<std::size_t>;


    struct DispatchTrack {

        /// The track parameter vector
        std::array<double, 6> parameters;

        /// The track covariance matrix
        std::array<double, 36> covariances;

        /// The indices of the measurements associated to this track
        std::vector<int> measurementIndices;

    };

} // namespace ActsExamples