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

// TODO make size configurable
std::array<float, 7 * 7> chargeMatrix(const ActsExamples::Cluster& cl) {
  // Compute the centroid index
  double acc0 = 0;
  double acc1 = 0;
  double tot = 0;
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    acc0 += cell.bin[0] * cell.activation;
    acc1 += cell.bin[1] * cell.activation;
    tot += cell.activation;
  }
  int icntr0 = static_cast<int>(acc0 / tot);
  int icntr1 = static_cast<int>(acc1 / tot);

  // By convention, put centroid in (3, 3) cell
  int diff0 = icntr0 - 3;
  int diff1 = icntr1 - 3;

  // Fill the matrix
  std::array<float, 7 * 7> mat{};
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    int im = cell.bin[0] - diff0;
    int jm = cell.bin[1] - diff1;
    if (im >= 0 && im < 7 && jm >= 0 && jm < 7)
      mat.at(im * 7 + jm) = cell.activation;
  }

  return mat;
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

  // TODO: Matrix size should be configurable perhaps?
  std::array<float, 7 * 7> matrix =
      ::detail::chargeMatrix((*clusters)[sourceLink.index()]);
  std::vector<float> input;
  input.assign(matrix.begin(), matrix.end());
  input.push_back(sourceLink.geometryId().volume());
  input.push_back(sourceLink.geometryId().layer());
  input.push_back(trackState.parameters()[Acts::eBoundPhi]);
  input.push_back(trackState.parameters()[Acts::eBoundTheta]);

  std::visit(
      [&](const auto& meas) {
        auto E = meas.expander();
        auto P = meas.projector();
        Acts::ActsVector<Acts::eBoundSize> fpar = E * meas.parameters();
        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * meas.covariance() * E.transpose();

        input.push_back(fpar[Acts::eBoundLoc0]);
        input.push_back(fpar[Acts::eBoundLoc1]);
        input.push_back(fcov(Acts::eBoundLoc0, Acts::eBoundLoc0));
        input.push_back(fcov(Acts::eBoundLoc1, Acts::eBoundLoc1));
        if (input.size() != n_inputs) {
          throw std::runtime_error("Expected input size of " +
                                   std::to_string(n_inputs) +
                                   ", got: " + std::to_string(input.size()));
        }

        std::vector<float> output = m_model.runONNXInference(input);
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
