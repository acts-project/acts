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
size_t fillChargeMatrix(Array& arr, const ActsExamples::Cluster& cl,
                        size_t size0 = 7u, size_t size1 = 7u) {
  // Compute the centroid index
  double acc0 = 0;
  double acc1 = 0;
  double tot = 0;
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    acc0 += cell.bin[0] * cell.activation;
    acc1 += cell.bin[1] * cell.activation;
    tot += cell.activation;
  }
  int icntr0 = (int)(acc0 / tot);
  int icntr1 = (int)(acc1 / tot);

  // By convention, put centroid in center of matrix
  // e.g. 3,3 for default matrix size of 7x7
  int diff0 = icntr0 - size0 / 2;
  int diff1 = icntr1 - size1 / 2;

  // Zero the charge matrix first, to guard against leftovers
  arr = Eigen::ArrayXXf::Zero(1, size0 * size1);
  // Fill the matrix
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    int im = cell.bin[0] - diff0;
    int jm = cell.bin[1] - diff1;
    if (im >= 0 && im < (int)size0 && jm >= 0 && jm < (int)size1) {
      typename Array::Index idx = im * size0 + jm;
      if (idx < arr.size()) {
        arr(idx) = cell.activation;
      }
    }
  }
  return size0 * size1;
}

}  // namespace detail

ActsExamples::NeuralCalibrator::NeuralCalibrator(
    const std::filesystem::path& modelPath, size_t nComp,
    std::vector<size_t> volumeIds)
    : m_env(ORT_LOGGING_LEVEL_WARNING, "NeuralCalibrator"),
      m_model(m_env, modelPath.c_str()),
      m_nComp{nComp},
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

  Acts::NetworkBatchInput inputBatch(1, n_inputs);
  auto input = inputBatch(0, Eigen::all);

  // TODO: Matrix size should be configurable perhaps?
  size_t matSize0 = 7u;
  size_t matSize1 = 7u;
  size_t i_input = ::detail::fillChargeMatrix(
      input, (*clusters)[sourceLink.index()], matSize0, matSize1);

  input[i_input++] = sourceLink.geometryId().volume();
  input[i_input++] = sourceLink.geometryId().layer();
  input[i_input++] = trackState.parameters()[Acts::eBoundPhi];
  input[i_input++] = trackState.parameters()[Acts::eBoundTheta];

  std::visit(
      [&](const auto& meas) {
        auto E = meas.expander();
        auto P = meas.projector();
        Acts::ActsVector<Acts::eBoundSize> fpar = E * meas.parameters();
        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * meas.covariance() * E.transpose();

        input[i_input++] = fpar[Acts::eBoundLoc0];
        input[i_input++] = fpar[Acts::eBoundLoc1];
        input[i_input++] = fcov(Acts::eBoundLoc0, Acts::eBoundLoc0);
        input[i_input++] = fcov(Acts::eBoundLoc1, Acts::eBoundLoc1);
        if (i_input != n_inputs) {
          throw std::runtime_error("Expected input size of " +
                                   std::to_string(n_inputs) +
                                   ", got: " + std::to_string(i_input));
        }

        // Input is a single row, hence .front()
        std::vector<float> output =
            m_model.runONNXInference(inputBatch).front();
        // Assuming 2-D measurements, the expected params structure is:
        // [      0,    nComp[ --> priors
        // [  nComp,  3*nComp[ --> means
        // [3*nComp,  5*nComp[ --> variances
        size_t nParams = 5 * m_nComp;
        if (output.size() != nParams) {
          throw std::runtime_error(
              "Got output vector of size " + std::to_string(output.size()) +
              ", expected size " + std::to_string(nParams));
        }

        // Most probable value computation of mixture density
        size_t iMax = 0;
        if (m_nComp > 1) {
          iMax = std::distance(
              output.begin(),
              std::max_element(output.begin(), output.begin() + m_nComp));
        }
        size_t iMpvLoc0 = m_nComp + iMax * 2;
        size_t iMpvVar0 = 3 * m_nComp + iMax * 2;

        fpar[Acts::eBoundLoc0] = output[iMpvLoc0];
        fpar[Acts::eBoundLoc1] = output[iMpvLoc0 + 1];
        fcov(Acts::eBoundLoc0, Acts::eBoundLoc0) = output[iMpvVar0];
        fcov(Acts::eBoundLoc1, Acts::eBoundLoc1) = output[iMpvVar0 + 1];

        constexpr size_t kSize =
            std::remove_reference_t<decltype(meas)>::size();
        std::array<Acts::BoundIndices, kSize> indices = meas.indices();
        Acts::ActsVector<kSize> cpar = P * fpar;
        Acts::ActsSymMatrix<kSize> ccov = P * fcov * P.transpose();

        Acts::SourceLink sl{sourceLink.geometryId(), sourceLink};

        Acts::Measurement<Acts::BoundIndices, kSize> cmeas(std::move(sl),
                                                           indices, cpar, ccov);

        trackState.allocateCalibrated(cmeas.size());
        trackState.setCalibrated(cmeas);
      },
      (measurements)[sourceLink.index()]);
}
