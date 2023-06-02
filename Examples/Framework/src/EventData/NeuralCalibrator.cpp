// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <ActsExamples/EventData/NeuralCalibrator.hpp>

#include <TFile.h>

namespace detail {

template <typename Array>
size_t fillChargeMatrix(Array& arr, const ActsExamples::Cluster& cluster,
                        size_t size0 = 7u, size_t size1 = 7u) {
  // Compute the centroid index
  double acc0 = 0;
  double acc1 = 0;
  double norm = 0;
  for (const ActsExamples::Cluster::Cell& cell : cluster.channels) {
    acc0 += cell.bin[0] * cell.activation;
    acc1 += cell.bin[1] * cell.activation;
    norm += cell.activation;
  }
  int iCenter0 = (int)(acc0 / norm);
  int iCenter1 = (int)(acc1 / norm);

  // By convention, put centroid in center of matrix
  // e.g. 3,3 for default matrix size of 7x7
  int offset0 = iCenter0 - size0 / 2;
  int offset1 = iCenter1 - size1 / 2;

  // Zero the charge matrix first, to guard against leftovers
  arr = Eigen::ArrayXXf::Zero(1, size0 * size1);
  // Fill the matrix
  for (const ActsExamples::Cluster::Cell& cell : cluster.channels) {
    int iMat = cell.bin[0] - offset0;
    int jMat = cell.bin[1] - offset1;
    if (iMat >= 0 && iMat < (int)size0 && jMat >= 0 && jMat < (int)size1) {
      typename Array::Index index = iMat * size0 + jMat;
      if (index < arr.size()) {
        arr(index) = cell.activation;
      }
    }
  }
  return size0 * size1;
}

}  // namespace detail

ActsExamples::NeuralCalibrator::NeuralCalibrator(
    const std::filesystem::path& modelPath, size_t nComponents,
    std::vector<size_t> volumeIds)
    : m_env(ORT_LOGGING_LEVEL_WARNING, "NeuralCalibrator"),
      m_model(m_env, modelPath.c_str()),
      m_nComponents{nComponents},
      m_volumeIds{std::move(volumeIds)} {}

void ActsExamples::NeuralCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& gctx,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
        trackState) const {
  Acts::SourceLink usl = trackState.getUncalibratedSourceLink();
  const IndexSourceLink& sourceLink = usl.get<IndexSourceLink>();
  assert((sourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  if (std::find(m_volumeIds.begin(), m_volumeIds.end(),
                sourceLink.geometryId().volume()) == m_volumeIds.end()) {
    m_fallback.calibrate(measurements, clusters, gctx, trackState);
    return;
  }

  Acts::NetworkBatchInput inputBatch(1, m_nInputs);
  auto input = inputBatch(0, Eigen::all);

  // TODO: Matrix size should be configurable perhaps?
  size_t matSize0 = 7u;
  size_t matSize1 = 7u;
  size_t iInput = ::detail::fillChargeMatrix(
      input, (*clusters)[sourceLink.index()], matSize0, matSize1);

  input[iInput++] = sourceLink.geometryId().volume();
  input[iInput++] = sourceLink.geometryId().layer();
  input[iInput++] = trackState.parameters()[Acts::eBoundPhi];
  input[iInput++] = trackState.parameters()[Acts::eBoundTheta];

  std::visit(
      [&](const auto& measurement) {
        auto E = measurement.expander();
        auto P = measurement.projector();
        Acts::ActsVector<Acts::eBoundSize> fpar = E * measurement.parameters();
        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * measurement.covariance() * E.transpose();

        input[iInput++] = fpar[Acts::eBoundLoc0];
        input[iInput++] = fpar[Acts::eBoundLoc1];
        input[iInput++] = fcov(Acts::eBoundLoc0, Acts::eBoundLoc0);
        input[iInput++] = fcov(Acts::eBoundLoc1, Acts::eBoundLoc1);
        if (iInput != m_nInputs) {
          throw std::runtime_error("Expected input size of " +
                                   std::to_string(m_nInputs) +
                                   ", got: " + std::to_string(iInput));
        }

        // Input is a single row, hence .front()
        std::vector<float> output =
            m_model.runONNXInference(inputBatch).front();
        // Assuming 2-D measurements, the expected params structure is:
        // [           0,    nComponent[ --> priors
        // [  nComponent,  3*nComponent[ --> means
        // [3*nComponent,  5*nComponent[ --> variances
        size_t nParams = 5 * m_nComponents;
        if (output.size() != nParams) {
          throw std::runtime_error(
              "Got output vector of size " + std::to_string(output.size()) +
              ", expected size " + std::to_string(nParams));
        }

        // Most probable value computation of mixture density
        size_t iMax = 0;
        if (m_nComponents > 1) {
          iMax = std::distance(
              output.begin(),
              std::max_element(output.begin(), output.begin() + m_nComponents));
        }
        size_t iLoc0 = m_nComponents + iMax * 2;
        size_t iVar0 = 3 * m_nComponents + iMax * 2;

        fpar[Acts::eBoundLoc0] = output[iLoc0];
        fpar[Acts::eBoundLoc1] = output[iLoc0 + 1];
        fcov(Acts::eBoundLoc0, Acts::eBoundLoc0) = output[iVar0];
        fcov(Acts::eBoundLoc1, Acts::eBoundLoc1) = output[iVar0 + 1];

        constexpr size_t kSize =
            std::remove_reference_t<decltype(measurement)>::size();
        std::array<Acts::BoundIndices, kSize> indices = measurement.indices();
        Acts::ActsVector<kSize> cpar = P * fpar;
        Acts::ActsSymMatrix<kSize> ccov = P * fcov * P.transpose();

        Acts::SourceLink sl{sourceLink.geometryId(), sourceLink};

        Acts::Measurement<Acts::BoundIndices, kSize> calibrated(
            std::move(sl), indices, cpar, ccov);

        trackState.allocateCalibrated(calibrated.size());
        trackState.setCalibrated(calibrated);
      },
      (measurements)[sourceLink.index()]);
}
