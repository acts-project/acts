// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/HistogramJsonConverter.hpp"

#include <boost/histogram/accumulators/mean.hpp>
#include <boost/histogram/algorithm/sum.hpp>
#include <boost/histogram/indexed.hpp>

using namespace Acts::Experimental;

namespace {

/// Serialize a single axis to JSON as {"edges":[...], "label":"..."}.
nlohmann::json axisToJson(const AxisVariant& axis) {
  std::vector<double> edges = extractBinEdges(axis);
  return nlohmann::json{{"edges", edges}, {"label", axis.metadata()}};
}

/// Build a JSON array of all axes from a boost histogram.
template <typename BH>
nlohmann::json axesJson(const BH& bh) {
  nlohmann::json axes = nlohmann::json::array();
  for (int i = 0; i < static_cast<int>(bh.rank()); ++i) {
    axes.push_back(axisToJson(bh.axis(i)));
  }
  return axes;
}

}  // namespace

nlohmann::json ActsPlugins::toJson(const Histogram1& boostHist) {
  const auto& bh = boostHist.histogram();

  std::vector<double> values;
  values.reserve(static_cast<std::size_t>(bh.axis(0).size()));
  for (auto&& x : boost::histogram::indexed(bh)) {
    values.push_back(static_cast<double>(*x));
  }

  return nlohmann::json{{"name", boostHist.name()},
                        {"title", boostHist.title()},
                        {"type", "histogram"},
                        {"axes", axesJson(bh)},
                        {"values", values}};
}

nlohmann::json ActsPlugins::toJson(const Histogram2& boostHist) {
  const auto& bh = boostHist.histogram();

  const int nx = bh.axis(0).size();
  const int ny = bh.axis(1).size();
  std::vector<double> values;
  values.resize(static_cast<std::size_t>(nx * ny), 0.0);
  for (auto&& x : boost::histogram::indexed(bh)) {
    values[static_cast<std::size_t>(x.index(0) * ny + x.index(1))] =
        static_cast<double>(*x);
  }

  return nlohmann::json{{"name", boostHist.name()},
                        {"title", boostHist.title()},
                        {"type", "histogram"},
                        {"axes", axesJson(bh)},
                        {"values", values}};
}

nlohmann::json ActsPlugins::toJson(const Histogram3& boostHist) {
  const auto& bh = boostHist.histogram();

  const int nx = bh.axis(0).size();
  const int ny = bh.axis(1).size();
  const int nz = bh.axis(2).size();
  std::vector<double> values;
  values.resize(static_cast<std::size_t>(nx * ny * nz), 0.0);
  for (auto&& x : boost::histogram::indexed(bh)) {
    values[static_cast<std::size_t>(x.index(0) * ny * nz + x.index(1) * nz +
                                    x.index(2))] = static_cast<double>(*x);
  }

  return nlohmann::json{{"name", boostHist.name()},
                        {"title", boostHist.title()},
                        {"type", "histogram"},
                        {"axes", axesJson(bh)},
                        {"values", values}};
}

nlohmann::json ActsPlugins::toJson(const ProfileHistogram1& boostProfile) {
  const auto& bh = boostProfile.histogram();
  const int n = bh.axis(0).size();

  using Accumulator = boost::histogram::accumulators::mean<double>;

  std::vector<double> counts, means, sumOfDeltasSquared;
  counts.reserve(static_cast<std::size_t>(n));
  means.reserve(static_cast<std::size_t>(n));
  sumOfDeltasSquared.reserve(static_cast<std::size_t>(n));

  for (auto&& x : boost::histogram::indexed(bh)) {
    const Accumulator& acc = *x;
    double count = acc.count();
    counts.push_back(count);
    means.push_back(count > 0 ? acc.value() : 0.0);
    sumOfDeltasSquared.push_back(count > 1.0 ? acc.variance() * (count - 1.0)
                                             : 0.0);
  }

  return nlohmann::json{{"name", boostProfile.name()},
                        {"title", boostProfile.title()},
                        {"type", "profile"},
                        {"sampleAxisTitle", boostProfile.sampleAxisTitle()},
                        {"axes", axesJson(bh)},
                        {"counts", counts},
                        {"means", means},
                        {"sum_of_deltas_squared", sumOfDeltasSquared}};
}

nlohmann::json ActsPlugins::toJson(const Efficiency1& boostEff) {
  const auto& accepted = boostEff.acceptedHistogram();
  const auto& total = boostEff.totalHistogram();
  const int n = accepted.axis(0).size();

  std::vector<double> acceptedVec, totalVec;
  acceptedVec.reserve(static_cast<std::size_t>(n));
  totalVec.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    acceptedVec.push_back(static_cast<double>(accepted.at(i)));
    totalVec.push_back(static_cast<double>(total.at(i)));
  }

  return nlohmann::json{{"name", boostEff.name()}, {"title", boostEff.title()},
                        {"type", "efficiency"},    {"axes", axesJson(accepted)},
                        {"accepted", acceptedVec}, {"total", totalVec}};
}

nlohmann::json ActsPlugins::toJson(const Efficiency2& boostEff) {
  const auto& accepted = boostEff.acceptedHistogram();
  const auto& total = boostEff.totalHistogram();
  const int nx = accepted.axis(0).size();
  const int ny = accepted.axis(1).size();

  std::vector<double> acceptedVec, totalVec;
  acceptedVec.resize(static_cast<std::size_t>(nx * ny), 0.0);
  totalVec.resize(static_cast<std::size_t>(nx * ny), 0.0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      std::size_t idx = static_cast<std::size_t>(i * ny + j);
      acceptedVec[idx] = static_cast<double>(accepted.at(i, j));
      totalVec[idx] = static_cast<double>(total.at(i, j));
    }
  }

  return nlohmann::json{{"name", boostEff.name()}, {"title", boostEff.title()},
                        {"type", "efficiency"},    {"axes", axesJson(accepted)},
                        {"accepted", acceptedVec}, {"total", totalVec}};
}
