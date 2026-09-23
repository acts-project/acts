// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"

#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Tolerance.hpp"

#include <format>

namespace {
    using namespace Acts::UnitLiterals;
    constexpr double inDeg(double angle) {
        return angle / 1._degree;
    }
}

namespace Acts::Experimental::detail {

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternState<Hit_t, Sector_t, Topology_t>::moveLineAnchorHit(
    const OrderedHit& refHit) 
{
    // Treat first the special case where we have only one group
    if (countGroups(/*onlyGoodGroups=*/ false) < 2) {
        // If we call this method with only one group, it means that we inverted the hit search direction without
        // finding any hit in other groups beside the initial one. So the anchor is the last added hit.
        lineAnchorHit = lastInsertedHit;
        return;
    }
    // Find first the closest group to the reference group among the pattern groups
    const auto& closestGroupIt = std::ranges::min_element(hitsPerGroup, std::ranges::less{},
        [&refHit](const auto& hits){
            if (hits.empty() || 
                    Topology_t::groupIndex(*hits.front()) == Topology_t::groupIndex(*refHit)) {
                return std::numeric_limits<int>::max();
            }
            return std::abs(hits.front().globLayer - refHit.globLayer);
        });

    // Then find the closest hit in that group to the reference hit
    const auto& hits {*closestGroupIt};
    lineAnchorHit = *std::ranges::min_element(hits, std::ranges::less{},
        [&refHit](const OrderedHit& hit){
            return std::abs(hit.globLayer - refHit.globLayer); });
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternState<Hit_t, Sector_t, Topology_t>::updateLineParameters(
    const GeometryContext& gctx,
    const BeamspotInfo& beamSpot) 
{
    if (!needLineUpdate) {
        return;
    }
    Vector3 pos1 {projToPhiPlane(gctx, *lineAnchorHit)};
    Vector3 pos2 {projToPhiPlane(gctx, *lastInsertedHit)};
    Vector3 d {pos2 - pos1};
    leverArm = d.norm();
    
    /** Check whether we have to use the beamspot instead of the anchor hit to draw the line. */
    useBeamspot = leverArm < cfg->minHitDistance4Line;
    if (useBeamspot) {
        linePos = beamSpot.position;
        d = pos2 - beamSpot.position;
        leverArm = d.norm();
    } else {
        linePos = pos1;
    }
    lineDir = d / leverArm;
    needLineUpdate = false;

    ACTS_VERBOSE(__func__<<"() Updated --> linePos: "<<toString(linePos)
        <<", lineDir: "<<toString(lineDir)<<", LeverArm: "
        <<leverArm<<", Use beamspot: "<<useBeamspot);
}
template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternState<Hit_t, Sector_t, Topology_t>::LineTestRes 
PatternState<Hit_t, Sector_t, Topology_t>::computeLineResidual(
    const GeometryContext& gctx, 
    const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const 
{
    LineTestRes res{};

    /** We project the test hit onto the phi plane only when the test hit does 
        *  not measure phi or when we have no phi layers, otherwise we do not project
        *  so the residual include the error in the phi direction. */
    const bool projectTestHit {!testHit->spacePoint()->measuresLoc0() || 
                               nPhiLayers == 0u};
    const Vector3 testPos {projectTestHit ? projToPhiPlane(gctx,*testHit) 
                                          : testHit->globalPosition(gctx)};
    const Vector3 K {testPos - linePos};
    const double KdotD {K.dot(lineDir)};
    
    const Vector3 R {K - KdotD * lineDir}; // Residual vector
    res.residual = R.norm();
    if (res.residual < Acts::s_epsilon) {
        /** If the residual is very small, it's likely due to bad topology, reject it. */
        res.residual = std::numeric_limits<double>::max();
        res.sigma = 0.;
        return res;
    }
    const Vector3 resDir {R / res.residual}; // Residual direction
    /** Alpha represents the extrapolation distance along the pattern line */
    const double alpha {KdotD / leverArm};
    
    /** Accumulate the derivatives of the scalar residual with respect to the common phi-plane angle. */
    double phiPlaneDerivativeAcc {0.};
    /** Accumulate the covariance contributions of the hits to the residual */
    double residualCovAcc {0.};

    /** @brief Accumulate the covariance contributions of one hit to the residual.
        *         Each hit contributes with its intrinsic covariance and, if projected,  
        *         with the uncertainty in the common phi-plane angle of the pattern.
        *   
        *  Intrinsic covariance: for a projected hit, the residual direction is  
        *  transformed with the projection jacobian J^T * dir.
        *  dir^T * J * cov * J^T * dir = (J^T * dir)^T * cov * (J^T * dir)
        *
        *  Phi plane uncertainty: for a projected hit, the derivative of the scalar  
        *  residual with respect to the common phi-plane angle for the projected hits.
        *
        *  @param hit: the hit for which to compute the covariance
        *  @param isProjected: flag indicating if the hit is projected
        *  @param pos: projected position (needed for the phi-plane derivative of projected hits)
        *  @param preFactor: factor, function of alpha, to scale the covariance contributions */
    auto covarianceTerm = [&](const Hit_t& hit,
                              const Vector3& pos,
                              const double preFactor,
                              bool isProjected) -> void {
                
        if (!isProjected) {
            residualCovAcc += Acts::square(preFactor) * 
                hit.intrinsicVariance(gctx, resDir);
            return;
        }
        const Vector3 sensorDir {hit.globalSensorDirection(gctx)};
        const double projFactor {sensorDir.dot(resDir) / 
                                    sensorDir.dot(bendPlaneNorm)};
        const Vector3 trfDir {resDir - projFactor * bendPlaneNorm};
        
        residualCovAcc += Acts::square(preFactor) 
            * hit.intrinsicVariance(gctx, trfDir);  
        phiPlaneDerivativeAcc += preFactor * Acts::fastHypot(pos.x(), pos.y()) * projFactor;
    };

    /** Compute the covariance contributions of the first line point */
    if (useBeamspot) {
        const double covS1 = beamSpot.length * Acts::square(resDir.z()) + 
                             beamSpot.radius * (1 - Acts::square(resDir.z()));
        residualCovAcc += Acts::square(alpha - 1) * covS1;
    } else {
        covarianceTerm(*lineAnchorHit, linePos, alpha - 1, /*isProjected=*/true);
    }

    /** Compute the covariance contributions of the second line point */
    const Vector3 pos2 {linePos + leverArm * lineDir};
    covarianceTerm(*lastInsertedHit, pos2, -alpha, /*isProjected=*/true);

    /** Compute the covariance contributions of the test hit */
    covarianceTerm(*testHit, testPos, 1., projectTestHit);
    
    res.sigma = std::sqrt(residualCovAcc + Acts::square(phiPlaneDerivativeAcc) * patPhiCov);

    ACTS_VERBOSE(__func__<<"() "<< brief(*this)<<"\nUse beamspot: "<<useBeamspot
        <<", alpha: "<<alpha<<", Residual: "<<res.residual<<" +- "<<res.sigma
        <<", linePos: "<<toString(linePos)<<", lineDir: "<<toString(lineDir)
        <<", testPos: "<<toString(testPos)<<", resDir: "<<toString(resDir)
        <<", hit pos sigma: "<<std::sqrt(residualCovAcc)
        <<", phi plane sigma: "<<std::abs(phiPlaneDerivativeAcc)*std::sqrt(patPhiCov));
    return res;
}
template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
Vector3
PatternState<Hit_t, Sector_t, Topology_t>::projToPhiPlane(
    const Acts::GeometryContext& gctx, 
    const Hit_t& hit) const 
{
    return Acts::PlanarHelper::intersectPlane(hit.globalPosition(gctx), 
                                              hit.globalSensorDirection(gctx),
                                              bendPlaneNorm, Vector3::Zero()).position();        
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternState<Hit_t, Sector_t, Topology_t>::updatePatternPhi(
    const GeometryContext& gctx)
{
    if (nPhiLayers == 0) {
        /** If there are no phi hits, we just use the central phi of the sector/overlap region, 
            *  with a standard deviation based on the expanded sector size. */
        patPhi = expSect.phi();
        patPhiCov = Acts::square(expSect.sectorSize()) / 3.;
        using namespace Acts::UnitLiterals;
        bendPlaneNorm = Acts::makeDirectionFromPhiTheta(patPhi + 90._degree, 90._degree);
        ACTS_VERBOSE(__func__<<"() No phi hits in the pattern, set pattern phi to "
            <<inDeg(patPhi)<<" +- "<<inDeg(std::sqrt(patPhiCov)));
        return;
    }
    double sumSin{0.}, sumCos{0.}, sumWeight{0.};

    auto processPhiHit = [&](const Hit_t& hit){
        if (!hit.spacePoint()->measuresLoc0()) {
            return;
        }
        const double phiCov {hit.phiVariance(gctx)};
        if (phiCov < Acts::s_epsilon) {
            std::stringstream ss {};
            ss << "Unexpected to have a phi hit with zero variance in phi direction: " << *hit.spacePoint() << "\n";
            throw std::runtime_error(ss.str());
        }
        const double w = 1./phiCov;

        const double phi {VectorHelpers::phi(hit.globalPosition(gctx))};
        sumSin += w * std::sin(phi);
        sumCos += w * std::cos(phi);
        sumWeight += w;
    };
    for (const std::vector<OrderedHit>&  hits : hitsPerGroup) {
        for (const auto& hit : hits) {
            processPhiHit(*hit);
        }
    }
    for (const Hit_t& hit : phiOnlyHits) {
        processPhiHit(hit);
    }
    
    patPhi = std::atan2(sumSin, sumCos);
    patPhiCov = 1./sumWeight;
    bendPlaneNorm = Acts::makeDirectionFromPhiTheta(patPhi + 90._degree, 90._degree);
    ACTS_VERBOSE(__func__<<"() Updated pattern phi to "
        <<inDeg(patPhi)<<" +- "<<inDeg(std::sqrt(patPhiCov)));
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
bool 
PatternState<Hit_t, Sector_t, Topology_t>::isPhiCompatible(
    const GeometryContext& gctx,
    const Hit_t& hit) const 
{
    /** We check that the test hit is compatible with the pattern phi, if available, which is given by the first
        *  phi measurement in the pattern. If the pattern doesn't have a phi yet, we check that the test hit is in 
        *  the same pattern sector(s) */
    const double testPhi {VectorHelpers::phi(hit.globalPosition(gctx))};
    if (nPhiLayers > 0) {
        const double deltaPhiSigma {std::sqrt(patPhiCov + hit.phiVariance(gctx))};
        const double deltaPhi {Acts::detail::difference_periodic(
            patPhi, testPhi, 2. * std::numbers::pi)};
        if (std::abs(deltaPhi) > cfg->nPhiSigma * deltaPhiSigma) {
            ACTS_VERBOSE(__func__<<"() The pattern with phi = "<<inDeg(patPhi)<<" +- "<<inDeg(std::sqrt(patPhiCov))
                <<" is not compatible with the test hit with phi "<<inDeg(testPhi) <<" +- "<<inDeg(std::sqrt(hit.phiVariance(gctx))));
            return false;
        }
    } else {
        if (!expSect.insideSector(testPhi)) {
            ACTS_VERBOSE(__func__<<"() The test hit with phi = "<<inDeg(testPhi)
                <<" is not inside the pattern sector: "<<static_cast<int>(expSect.sector()));
            return false;
        }            
    }
    return true;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternState<Hit_t, Sector_t, Topology_t>::addHit(
    const GeometryContext& gctx,
    const OrderedHit& hit,
    const double residual,
    const double resSigma) 
{
    /** Add the new hit */
    hitsPerGroup[Topology_t::groupIndex(*hit)].push_back(hit);

    /** Update the hit counts */
    if (hit->isPrecision()) {
        nPrecisionLayers++;
    } else {
        nTriggerLayers++;
    }
    if (hit->spacePoint()->measuresLoc0()) {
        nPhiLayers++;
        updatePatternPhi(gctx);
    }

    /** Update the pointers to previous layer hit */
    const bool isNewSGroup {Topology_t::groupIndex(*hit) != 
                            Topology_t::groupIndex(*lastInsertedHit)};
    prevLayerHit = lastInsertedHit;
    lastInsertedHit = hit;

    meanNormResidual2 += Acts::square(residual / resSigma);
    lastResSigma = resSigma;
    lastResidual = residual;

    /** If the new compatible hit is in a different group, update the line anchor */
    if (isNewSGroup) {
        moveLineAnchorHit(hit);
    }
    needLineUpdate = true;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternState<Hit_t, Sector_t, Topology_t>::overWriteHit(
    const GeometryContext& gctx,
    const OrderedHit& newHit,
    const double newResidual,
    const double newResSigma) 
{
    const auto group {Topology_t::groupIndex(*newHit)};
    if (group != Topology_t::groupIndex(*lastInsertedHit) || 
            lastInsertedHit.globLayer != newHit.globLayer) {
        throw std::runtime_error(std::format(
            "Trying to overwrite a hit in group/layer {}/{} with another one from group/layer {}/{}", 
                static_cast<int>(Topology_t::groupIndex(*lastInsertedHit)), 
                lastInsertedHit.globLayer, static_cast<int>(group), newHit.globLayer));
    }
    /* We expect to overwrite hits of the same type (precision/trigger), since we only branch when we have 
        * compatible hits in the same layer, except for sTGC hits, where we have pad and strips in the same layer */
    if (lastInsertedHit->isPrecision() != newHit->isPrecision()) {
        if (newHit->isPrecision()) {
            nPrecisionLayers++;
            nTriggerLayers--;
        } else {
            nPrecisionLayers--;
            nTriggerLayers++;
        }
    }
    /** Update the phi counts */
    bool updatePhi {false};
    if (lastInsertedHit->spacePoint()->measuresLoc0()) {
        nPhiLayers--;
        updatePhi = true;
    }
    if (newHit->spacePoint()->measuresLoc0()) {
        nPhiLayers++;
        updatePhi = true;
    }
    /** Update the residual */
    meanNormResidual2 += Acts::square(newResidual / newResSigma) - 
                         Acts::square(lastResidual / lastResSigma);
    lastResSigma = newResSigma;
    lastResidual = newResidual;

    auto& stHits {hitsPerGroup[group]};
    if (stHits.back() != lastInsertedHit) {
        std::stringstream ss {};
        ss << "Trying to overwrite a hit that is not the last inserted hit in group/layer " 
        << static_cast<int>(group) << "/" << lastInsertedHit.globLayer << "\n";
        ss << "Last inserted hit: " << *lastInsertedHit->spacePoint() << "\n";
        ss << "Last hit in group: " << *stHits.back()->spacePoint();
        throw std::runtime_error(ss.str());
    }
    stHits.pop_back();

    /** Add the new hit */
    stHits.push_back(newHit);
    lastInsertedHit = newHit;

    if (updatePhi) {
        updatePatternPhi(gctx);
    }
    needLineUpdate = true;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
bool 
PatternState<Hit_t, Sector_t, Topology_t>::isInPattern(
    const Hit_t& hit) const 
{
    const auto& hits {hitsPerGroup[Topology_t::groupIndex(hit)]};
    return std::ranges::find_if(hits, 
        [&hit](const OrderedHit& c){ return *c == hit; }) != hits.end();                                           
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
double 
PatternState<Hit_t, Sector_t, Topology_t>::getMeanResidual2() const {
    if (isFinalized) {
        return meanNormResidual2;
    }
    return meanNormResidual2 / nBendingLayers();
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename Topology_t::LayerIdx 
PatternState<Hit_t, Sector_t, Topology_t>::nBendingLayers() const {
    return nPrecisionLayers + nTriggerLayers;
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename Topology_t::GroupIdx 
PatternState<Hit_t, Sector_t, Topology_t>::countGroups(
    const bool onlyGoodGroups) const 
{
    typename Topology_t::GroupIdx nGroups = 0u;
    for (typename Topology_t::GroupIdx g = 0u; g < Topology_t::nGroups; ++g) {
        const auto& hits {hitsPerGroup[g]};
        if (!hits.empty() && 
                (!onlyGoodGroups || hits.size() >= cfg->minGroupLayers)) {
            nGroups++;
        }
    }
    return nGroups;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternState<Hit_t, Sector_t, Topology_t>::PatternPrintView 
PatternState<Hit_t, Sector_t, Topology_t>::brief(const PatternState& p) {
    return {p, /*detailed=*/false};
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternState<Hit_t, Sector_t, Topology_t>::PatternPrintView 
PatternState<Hit_t, Sector_t, Topology_t>::detailed(const PatternState& p) {
    return {p, /*detailed=*/true};
}


template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
PatternState<Hit_t, Sector_t, Topology_t>::PatternState(
    const GeometryContext& gctx,
    const OrderedHit& seed,
    const typename Sector_t::Index_t  expSector,          
    const Config* cfg,
    const Acts::Logger* logger)
        : cfg{cfg},
          m_logger{logger},
          lastInsertedHit{seed},
          prevLayerHit{seed},
          lineAnchorHit{seed},
          patTheta{VectorHelpers::theta(seed->globalPosition(gctx))},
          expSect{expSector} {
                
        /** Add the new hit */
        hitsPerGroup[Topology_t::groupIndex(*seed)].push_back(seed);

        /** Update the hit counts */
        if (seed->isPrecision()) {
            nPrecisionLayers++;
        } else {
            nTriggerLayers++;
        }
        if (seed->spacePoint()->measuresLoc0()) {
            nPhiLayers++;
        }
        updatePatternPhi(gctx);
        needLineUpdate = true;
    }

template<GlobPatFinderHit Hit_t, 
         SectorType Sector_t, 
         PatternTopology<Hit_t> Topology_t>
void PatternState<Hit_t, Sector_t, Topology_t>::print(
    std::ostream& ostr, 
    bool detailed) const 
{
    using namespace Acts::UnitLiterals;
    auto inDeg = [](double angle) { return angle / 1._degree; };
    ostr<<"PatternState Exp Sector: "<<static_cast<int>(expSect.sector())
    <<", Theta: "<<inDeg(patTheta) << ", Phi: "<<inDeg(patPhi)<<" +- "<<inDeg(std::sqrt(patPhiCov));
    ostr<<", nPrec: "<<static_cast<int>(nPrecisionLayers)<<", nEtaNonPrec: "
        <<static_cast<int>(nTriggerLayers)<<", nPhi: "<<static_cast<int>(nPhiLayers);
    ostr<<", mean norma res sq: "<<getMeanResidual2()<<", dir: "<<toString(lineDir);
    ostr<<", Hit per group: \n";
    for (typename Topology_t::GroupIdx g = 0u; g < Topology_t::nGroups; ++g) {
        const auto& hits {hitsPerGroup[g]};
        if (hits.empty()) {
            continue;
        }

        ostr<<"  Group "<<static_cast<int>(g)<<" has "<<hits.size()<<" hits ";
        if (detailed) {
            ostr<<"\n";
            for (const auto& hit : hits) {
                ostr<<"    "<<*hit->spacePoint()<<", globLay: "<<static_cast<int>(hit.globLayer)<<"\n";
            }
        }
    }
    if (!detailed) {
        ostr<<"\n";
        ostr<<"  Last inserted hit: "<<lastInsertedHit<<"\n";
        ostr<<"  Previous layer hit: "<<prevLayerHit<<"\n";
        ostr<<"  Line anchor hit: "<<lineAnchorHit<<"\n";
    }
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void
PatternState<Hit_t, Sector_t, Topology_t>::OrderedHit::print(std::ostream& ostr) const {
    ostr<<*hitPtr->spacePoint()<< ", group: " << static_cast<int>(Topology_t::groupIndex(*hitPtr)) 
        <<", globLay: "<<static_cast<int>(globLayer);
}
}