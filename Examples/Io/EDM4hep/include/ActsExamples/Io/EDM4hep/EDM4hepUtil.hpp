// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <functional>

#include "edm4hep/MCParticle.h"
#include "edm4hep/MutableMCParticle.h"
#include "edm4hep/MutableSimTrackerHit.h"
#include "edm4hep/MutableTrack.h"
#include "edm4hep/MutableTrackerHit.h"
#include "edm4hep/MutableTrackerHitPlane.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlane.h"

namespace ActsExamples {
namespace EDM4hepUtil {

static constexpr std::int32_t EDM4HEP_ACTS_POSITION_TYPE = 42;

using MapParticleIdFrom =
    std::function<ActsFatras::Barcode(edm4hep::MCParticle particle)>;
using MapParticleIdTo =
    std::function<edm4hep::MCParticle(ActsFatras::Barcode particleId)>;

using MapGeometryIdFrom =
    std::function<Acts::GeometryIdentifier(std::uint64_t cellId)>;
using MapGeometryIdTo =
    std::function<std::uint64_t(Acts::GeometryIdentifier geometryId)>;

/// Read a Fatras particle from EDM4hep.
///
/// Note: The particle ID is not persitent in EDM4hep.
ActsFatras::Particle readParticle(edm4hep::MCParticle from,
                                  MapParticleIdFrom particleMapper);

/// Write a Fatras particle into EDM4hep.
///
/// Note: The particle ID is not persitent in EDM4hep.
void writeParticle(const ActsFatras::Particle& from,
                   edm4hep::MutableMCParticle to);

/// Read a Fatras hit from EDM4hep.
///
/// Note: The hit index is not persitent.
ActsFatras::Hit readSimHit(const edm4hep::SimTrackerHit& from,
                           MapParticleIdFrom particleMapper,
                           MapGeometryIdFrom geometryMapper);

/// Write a Fatras hit to EDM4hep.
///
/// Note: The hit index is not persitent.
void writeSimHit(const ActsFatras::Hit& from, edm4hep::MutableSimTrackerHit to,
                 MapParticleIdTo particleMapper,
                 MapGeometryIdTo geometryMapper);

/// Read a measurement cluster from EDM4hep.
///
/// Note: The hit index is not persitent.
/// Note: 1D not supported.
/// Note: Clusters are written to inappropriate fields.
/// Note: Digitization is not supported.
Measurement readMeasurement(edm4hep::TrackerHitPlane from,
                            const edm4hep::TrackerHitCollection* fromClusters,
                            Cluster* toCluster,
                            MapGeometryIdFrom geometryMapper);

/// Write a measurement cluster to EDM4hep.
///
/// Note: The hit index is not persitent.
/// Note: 1D not supported.
/// Note: Clusters are written to inappropriate fields.
/// Note: Digitization is not supported.
void writeMeasurement(const Measurement& from,
                      edm4hep::MutableTrackerHitPlane to,
                      const Cluster* fromCluster,
                      edm4hep::TrackerHitCollection& toClusters,
                      MapGeometryIdTo geometryMapper);

void writeTrajectory(const Trajectories& from, edm4hep::MutableTrack to,
                     std::size_t fromIndex,
                     const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap);

}  // namespace EDM4hepUtil
}  // namespace ActsExamples
