// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <numbers>

#include <edm4hep/EDM4hepVersion.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/MutableVertex.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/Vector4f.h>

using namespace Acts;
using namespace Acts::detail;

namespace ActsPlugins::EDM4hepUtil {
namespace detail {

SquareMatrix<6> jacobianToEdm4hep(double theta, double qOverP, double Bz) {
  // Calculate jacobian from our internal parametrization (d0, z0, phi, theta,
  // q/p) to the LCIO / edm4hep (see:
  // https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf)
  // one (d0, z0, phi, tan(lambda), omega). Top left 3x3 matrix in the
  // jacobian is 1. Bottom right 2x2 matrix is:
  //
  // [  d                                 ]
  // [------(cot(theta))         0        ]
  // [dtheta                              ]
  // [                                    ]
  // [  d   /  B*q/p   \   d  /  B*q/p   \]
  // [------|----------|  ----|----------|]
  // [dtheta\sin(theta)/  dq/p\sin(theta)/]
  //
  // =
  //
  // [     2                        ]
  // [- csc (theta)           0     ]
  // [                              ]
  // [-B*q/p*cos(theta)       B     ]
  // [------------------  ----------]
  // [      2             sin(theta)]
  // [   sin (theta)                ]

  SquareMatrix<6> J;
  J.setIdentity();
  double sinTheta = std::sin(theta);
  J(3, 3) = -1.0 / (sinTheta * sinTheta);
  J(4, 4) = Bz / sinTheta;  // dOmega / d(qop)
  J(4, 3) = -Bz * qOverP * std::cos(theta) /
            (sinTheta * sinTheta);  // dOmega / dTheta
  return J;
}

SquareMatrix<6> jacobianFromEdm4hep(double tanLambda, double omega, double Bz) {
  // [     d      /                     pi\                                  ]
  // [------------|-atan(\tan\lambda) + --|                 0                ]
  // [d\tan\lambda\                     2 /                                  ]
  // [                                                                       ]
  // [     d      /         \Omega        \     d   /         \Omega        \]
  // [------------|-----------------------|  -------|-----------------------|]
  // [d\tan\lambda|     __________________|  d\Omega|     __________________|]
  // [            |    /            2     |         |    /            2     |]
  // [            \B*\/  \tan\lambda  + 1 /         \B*\/  \tan\lambda  + 1 /]
  //
  // =
  //
  // [         -1                                     ]
  // [   ----------------                 0           ]
  // [              2                                 ]
  // [   \tan\lambda  + 1                             ]
  // [                                                ]
  // [  -\Omega*\tan\lambda               1           ]
  // [-----------------------  -----------------------]
  // [                    3/2       __________________]
  // [  /           2    \         /            2     ]
  // [B*\\tan\lambda  + 1/     B*\/  \tan\lambda  + 1 ]

  SquareMatrix<6> J;
  J.setIdentity();
  J(3, 3) = -1 / (tanLambda * tanLambda + 1);
  J(4, 3) = -1 * omega * tanLambda /
            (Bz * std::pow(tanLambda * tanLambda + 1, 3. / 2.));
  J(4, 4) = 1 / (Bz * std::hypot(tanLambda, 1));

  return J;
}

void packCovariance(const SquareMatrix<6>& from, float* to) {
  for (int i = 0; i < from.rows(); i++) {
    for (int j = 0; j <= i; j++) {
      std::size_t k = (i + 1) * i / 2 + j;
      to[k] = static_cast<float>(from(i, j));
    }
  }
}

void unpackCovariance(const float* from, SquareMatrix<6>& to) {
  auto k = [](std::size_t i, std::size_t j) { return (i + 1) * i / 2 + j; };
  for (int i = 0; i < to.rows(); i++) {
    for (int j = 0; j < to.cols(); j++) {
      to(i, j) = from[j <= i ? k(i, j) : k(j, i)];
    }
  }
}

Parameters convertTrackParametersToEdm4hep(const GeometryContext& gctx,
                                           double Bz,
                                           const BoundTrackParameters& params) {
  Vector3 global = params.referenceSurface().localToGlobal(
      gctx, params.parameters().template head<2>(), params.direction());

  std::shared_ptr<const Surface> refSurface =
      params.referenceSurface().getSharedPtr();
  BoundVector targetPars = params.parameters();
  std::optional<BoundMatrix> targetCov = params.covariance();

  // If the reference surface is a perigee surface, we use that. Otherwise
  // we create a new perigee surface at the global position of the track
  // parameters.
  if (dynamic_cast<const PerigeeSurface*>(refSurface.get()) == nullptr) {
    refSurface = Surface::makeShared<PerigeeSurface>(global);

    // We need to convert to the target parameters
    // Keep the free parameters around we might need them for the covariance
    // conversion

    auto perigeeParams =
        boundToBoundConversion(gctx, params, *refSurface, Vector3{0, 0, Bz})
            .value();
    targetPars = perigeeParams.parameters();
    targetCov = perigeeParams.covariance();
  }

  Parameters result;
  result.surface = refSurface;

  // Only run covariance conversion if we have a covariance input
  if (targetCov) {
    SquareMatrix<6> J = jacobianToEdm4hep(targetPars[eBoundTheta],
                                          targetPars[eBoundQOverP], Bz);
    result.covariance = J * targetCov.value() * J.transpose();
  }

  result.values[0] = targetPars[eBoundLoc0];
  result.values[1] = targetPars[eBoundLoc1];
  result.values[2] = targetPars[eBoundPhi];
  result.values[3] = std::tan(std::numbers::pi / 2. - targetPars[eBoundTheta]);
  result.values[4] =
      targetPars[eBoundQOverP] / std::sin(targetPars[eBoundTheta]) * Bz;
  result.values[5] = targetPars[eBoundTime];

  result.particleHypothesis = params.particleHypothesis();

  return result;
}

BoundTrackParameters convertTrackParametersFromEdm4hep(
    double Bz, const Parameters& params) {
  BoundVector targetPars;

  SquareMatrix<6> J =
      jacobianFromEdm4hep(params.values[3], params.values[4], Bz);

  std::optional<BoundMatrix> cov;
  if (params.covariance.has_value()) {
    cov = J * params.covariance.value() * J.transpose();
  }

  targetPars[eBoundLoc0] = params.values[0];
  targetPars[eBoundLoc1] = params.values[1];
  targetPars[eBoundPhi] = params.values[2];
  targetPars[eBoundTheta] = std::numbers::pi / 2. - std::atan(params.values[3]);
  targetPars[eBoundQOverP] =
      params.values[4] * std::sin(targetPars[eBoundTheta]) / Bz;
  targetPars[eBoundTime] = params.values[5];

  return {params.surface, targetPars, cov, params.particleHypothesis};
}

}  // namespace detail

#if EDM4HEP_VERSION_MAJOR >= 1 || \
    (EDM4HEP_VERSION_MAJOR == 0 && EDM4HEP_VERSION_MINOR == 99)
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit) {
  return hit.getParticle();
}

void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle) {
  hit.setParticle(particle);
}
#else
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit) {
  return hit.getMCParticle();
}

void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle) {
  hit.setMCParticle(particle);
}
#endif

std::size_t SimHitAssociation::size() const {
  return m_internalToEdm4hep.size();
}

void SimHitAssociation::reserve(std::size_t size) {
  m_internalToEdm4hep.reserve(size);
}

void SimHitAssociation::add(std::size_t internalIndex,
                            const edm4hep::SimTrackerHit& edm4hepHit) {
  m_internalToEdm4hep.push_back(edm4hepHit);
  // m_edm4hepToInternal.at(edm4hepHit.id()) = internalIndex;
  m_edm4hepToInternal.emplace(edm4hepHit.id(), internalIndex);
}

edm4hep::SimTrackerHit SimHitAssociation::lookup(
    std::size_t internalIndex) const {
  return m_internalToEdm4hep.at(internalIndex);
}

std::size_t SimHitAssociation::lookup(const edm4hep::SimTrackerHit& hit) const {
  return m_edm4hepToInternal.at(hit.id());
}

void writeVertex(const Vertex& vertex, edm4hep::MutableVertex to) {
  static constexpr std::array<edm4hep::FourMomCoords, 4> toEdm4hep = []() {
    std::array<edm4hep::FourMomCoords, 4> values{};
    using enum edm4hep::FourMomCoords;
    values.at(eFreePos0) = x;
    values.at(eFreePos1) = y;
    values.at(eFreePos2) = z;
    values.at(eFreeTime) = t;
    return values;
  }();

  // Wrap this in a templated lambda so we can use `if constexpr` to select the
  // correct write function based on the type of the properties of the `to`
  // object.
  auto writeVertex = []<typename T>(const Vertex& vertex_, T& to_)
    requires(std::is_same_v<T, edm4hep::MutableVertex>)
  {
    if constexpr (detail::kEdm4hepVertexHasTime) {
      Vector4 pos = vertex_.fullPosition();
      to_.setPosition({static_cast<float>(pos[eFreePos0]),
                       static_cast<float>(pos[eFreePos1]),
                       static_cast<float>(pos[eFreePos2]),
                       static_cast<float>(pos[eFreeTime])});

      edm4hep::CovMatrix4f& cov = to_.getCovMatrix();
      std::array coords{eFreePos0, eFreePos1, eFreePos2, eFreeTime};
      for (auto i : coords) {
        for (auto j : coords) {
          cov.setValue(static_cast<float>(vertex_.fullCovariance()(i, j)),
                       toEdm4hep.at(i), toEdm4hep.at(j));
        }
      }
    } else {
      Vector3 pos = vertex_.position();
      to_.setPosition({static_cast<float>(pos[eFreePos0]),
                       static_cast<float>(pos[eFreePos1]),
                       static_cast<float>(pos[eFreePos2])});
      edm4hep::CovMatrix3f& cov = to_.getCovMatrix();
      std::array coords{eFreePos0, eFreePos1, eFreePos2};
      for (auto i : coords) {
        for (auto j : coords) {
          cov.setValue(static_cast<float>(vertex_.covariance()(i, j)),
                       toEdm4hep.at(i), toEdm4hep.at(j));
        }
      }
    }

    to_.setChi2(static_cast<float>(vertex_.fitQuality().first));
    to_.setNdf(static_cast<int>(vertex_.fitQuality().second));
  };

  writeVertex(vertex, to);
}

namespace detail {
/// Encode a list of bound parameter indices into a 32-bit integer.
///
/// This function bit-packs up to 6 parameter indices (each 0-5 corresponding to
/// eBoundLoc0, eBoundLoc1, eBoundPhi, eBoundTheta, eBoundQOverP, eBoundTime)
/// into a single 32-bit unsigned integer for storage in EDM4hep format.
///
/// Bit layout:
///   - Bits 0-3:   Number of indices (size)
///   - Bits 4-7:   First index
///   - Bits 8-11:  Second index
///   - Bits 12-15: Third index
///   - (and so on, up to 6 indices total)
///
/// @param indices Span of parameter indices to encode (max 6 elements)
/// @return Packed 32-bit unsigned integer containing all indices
std::uint32_t encodeIndices(std::span<const std::uint8_t> indices) {
  if (indices.size() > eBoundSize) {
    throw std::runtime_error(
        "Number of indices exceeds maximum of 6 for EDM4hep");
  }
  std::uint32_t result = 0;

  std::uint8_t shift = 0;
  result |= (indices.size() << 0);
  shift += 4;

  for (std::uint8_t index : indices) {
    if (index > eBoundSize) {
      throw std::runtime_error(
          "Index out of range: can only encode indices up to 4 bits (0-15)");
    }
    result |= (static_cast<std::uint32_t>(index) << shift);
    shift += 4;
  }
  return result;
}

/// Decode a 32-bit integer back into a list of bound parameter indices.
///
/// This function unpacks a bit-packed integer (created by encodeIndices) back
/// into the original list of parameter indices. See encodeIndices for the bit
/// layout specification.
///
/// @param type Packed 32-bit unsigned integer containing encoded indices
/// @return Vector of decoded parameter indices (0-5 for each bound parameter)
boost::container::static_vector<std::uint8_t, eBoundSize> decodeIndices(
    std::uint32_t type) {
  boost::container::static_vector<std::uint8_t, eBoundSize> result;
  std::uint8_t size = type & 0xF;
  if (size > eBoundSize) {
    throw std::runtime_error(
        "Number of indices exceeds maximum of 6 for EDM4hep");
  }
  result.resize(size);
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = (type >> ((i + 1) * 4)) & 0xF;
    if (result[i] > eBoundSize) {
      throw std::runtime_error(
          "Index out of range: can only encode indices up to 4 bits (0-15)");
    }
  }
  return result;
}
}  // namespace detail

void writeMeasurement(const GeometryContext& gctx,
                      const Eigen::Map<const ActsDynamicVector>& parameters,
                      const Eigen::Map<const ActsDynamicMatrix>& covariance,
                      std::span<const std::uint8_t> indices,
                      std::uint64_t cellId, const Acts::Surface& surface,
                      ActsPodioEdm::MutableTrackerHitLocal to) {
  if (parameters.size() != covariance.rows() ||
      covariance.rows() != covariance.cols() || parameters.size() < 0 ||
      indices.size() != static_cast<std::size_t>(parameters.size())) {
    throw std::runtime_error(
        "Size mismatch between parameters and covariance matrix");
  }

  auto dim = static_cast<std::size_t>(parameters.size());

  if (cellId != 0) {
    to.setCellID(cellId);
  }

  to.setType(detail::encodeIndices(indices));

  auto loc0 = std::ranges::find(indices, eBoundLoc0);
  auto loc1 = std::ranges::find(indices, eBoundLoc1);
  auto time = std::ranges::find(indices, eBoundTime);

  if (loc0 != indices.end() && loc1 != indices.end()) {
    Vector2 loc{parameters[std::distance(indices.begin(), loc0)],
                parameters[std::distance(indices.begin(), loc1)]};
    Vector3 global = surface.localToGlobal(gctx, loc, Vector3::UnitZ());
    global /= Acts::UnitConstants::mm;
    to.setPosition({global.x(), global.y(), global.z()});
  }

  if (time != indices.end()) {
    std::size_t timeOffset = std::distance(indices.begin(), time);
    to.setTime(
        static_cast<float>(parameters[timeOffset] / Acts::UnitConstants::ns));
  }

  for (double value : std::span{parameters.data(), dim}) {
    to.addToMeasurement(static_cast<float>(value));
  }

  for (double value : std::span{covariance.data(), dim * dim}) {
    to.addToCovariance(static_cast<float>(value));
  }
}

}  // namespace ActsPlugins::EDM4hepUtil
