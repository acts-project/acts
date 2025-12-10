// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackSelector.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;

namespace ActsPython {

/// Add the track finding related bindings to the module
/// @param m The module to add the bindings to
void addTrackFinding(py::module_& m) {
  {
    auto constructor =
        [](const std::vector<std::pair<
               Acts::GeometryIdentifier,
               std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, std::vector<std::size_t>>>>&
               input) {
          std::vector<std::pair<Acts::GeometryIdentifier,
                                Acts::MeasurementSelectorCuts>>
              converted;
          converted.reserve(input.size());
          for (const auto& [id, cuts] : input) {
            const auto& [bins, chi2Measurement, chi2Outlier, num] = cuts;
            converted.emplace_back(
                id, Acts::MeasurementSelectorCuts{bins, chi2Measurement, num,
                                                  chi2Outlier});
          }
          return std::make_unique<Acts::MeasurementSelector::Config>(converted);
        };

    py::class_<Acts::MeasurementSelectorCuts>(m, "MeasurementSelectorCuts")
        .def(py::init<>())
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<std::size_t>, std::vector<double>>())
        .def_readwrite("etaBins", &Acts::MeasurementSelectorCuts::etaBins)
        .def_readwrite("chi2CutOffMeasurement",
                       &Acts::MeasurementSelectorCuts::chi2CutOff)
        .def_readwrite("chi2CutOffOutlier",
                       &Acts::MeasurementSelectorCuts::chi2CutOffOutlier)
        .def_readwrite("numMeasurementsCutOff",
                       &Acts::MeasurementSelectorCuts::numMeasurementsCutOff);

    auto ms = py::class_<Acts::MeasurementSelector>(m, "MeasurementSelector");
    auto c = py::class_<Acts::MeasurementSelector::Config>(ms, "Config")
                 .def(py::init<
                      std::vector<std::pair<Acts::GeometryIdentifier,
                                            Acts::MeasurementSelectorCuts>>>())
                 .def(py::init(constructor));
  }
  {
    using EtaBinnedConfig = TrackSelector::EtaBinnedConfig;
    using Config = TrackSelector::Config;

    auto tool = py::class_<TrackSelector>(m, "TrackSelector")
                    .def(py::init<const Config&>(), py::arg("config"))
                    .def(py::init<const EtaBinnedConfig&>(), py::arg("config"));

    {
      auto mc = py::class_<TrackSelector::MeasurementCounter>(
                    tool, "MeasurementCounter")
                    .def(py::init<>())
                    .def("addCounter",
                         &TrackSelector::MeasurementCounter::addCounter);
    }

    {
      auto c = py::class_<Config>(tool, "Config").def(py::init<>());

      patchKwargsConstructor(c);

      ACTS_PYTHON_STRUCT(c, loc0Min, loc0Max, loc1Min, loc1Max, timeMin,
                         timeMax, phiMin, phiMax, etaMin, etaMax, absEtaMin,
                         absEtaMax, ptMin, ptMax, minMeasurements, maxHoles,
                         maxOutliers, maxHolesAndOutliers, maxSharedHits,
                         maxChi2, measurementCounter, requireReferenceSurface);

      pythonRangeProperty(c, "loc0", &Config::loc0Min, &Config::loc0Max);
      pythonRangeProperty(c, "loc1", &Config::loc1Min, &Config::loc1Max);
      pythonRangeProperty(c, "time", &Config::timeMin, &Config::timeMax);
      pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
      pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
      pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
      pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
    }

    {
      auto c = py::class_<EtaBinnedConfig>(tool, "EtaBinnedConfig")
                   .def(py::init<>())
                   .def(py::init<const Config&>());

      patchKwargsConstructor(c);

      c.def_property_readonly("nEtaBins", &EtaBinnedConfig::nEtaBins);

      ACTS_PYTHON_STRUCT(c, cutSets, absEtaEdges);
    }
  }
}
}  // namespace ActsPython
