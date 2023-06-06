// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Definitions/TrackParametrization.hpp>
#include <ActsExamples/EventData/ScalingCalibrator.hpp>

#include <regex>

#include <TKey.h>

namespace detail {

std::pair<Acts::GeometryIdentifier, std::string> parseMapKey(
    const std::string& mapkey) {
  std::regex reg("^map_([0-9]+)-([0-9]+)-([0-9]+)_([xy]_.*)$");
  std::smatch matches;

  if (std::regex_search(mapkey, matches, reg) && matches.size() == 5) {
    size_t vol = std::stoull(matches[1].str());
    size_t lyr = std::stoull(matches[2].str());
    size_t mod = std::stoull(matches[3].str());

    Acts::GeometryIdentifier geoId;
    geoId.setVolume(vol);
    geoId.setLayer(lyr);
    geoId.setSensitive(mod);

    std::string var(matches[4].str());

    return std::make_pair(geoId, var);
  } else {
    throw std::runtime_error("Invalid map key: " + mapkey);
  }
}

std::map<Acts::GeometryIdentifier, ActsExamples::ScalingCalibrator::MapTuple>
readMaps(const std::filesystem::path& path) {
  std::map<Acts::GeometryIdentifier, ActsExamples::ScalingCalibrator::MapTuple>
      maps;

  TFile ifile(path.c_str(), "READ");
  if (ifile.IsZombie()) {
    throw std::runtime_error("Unable to open TFile: " + path.string());
  }

  TList* lst = ifile.GetListOfKeys();
  assert(lst != nullptr);

  for (auto it = lst->begin(); it != lst->end(); ++it) {
    TKey* key = static_cast<TKey*>(*it);
    if (std::strcmp(key->GetClassName(), "TH2D") == 0) {
      auto [geoId, var] = parseMapKey(key->GetName());

      TH2D hist;
      key->Read(&hist);

      if (var == "x_offset") {
        maps[geoId].x_offset = hist;
      } else if (var == "x_scale") {
        maps[geoId].x_scale = hist;
      } else if (var == "y_offset") {
        maps[geoId].y_offset = hist;
      } else if (var == "y_scale") {
        maps[geoId].y_scale = hist;
      } else {
        throw std::runtime_error("Unrecognized var: " + var);
      }
    }
  }
  return maps;
}

std::bitset<3> readMask(const std::filesystem::path& path) {
  TFile ifile(path.c_str(), "READ");
  if (ifile.IsZombie()) {
    throw std::runtime_error("Unable to open TFile: " + path.string());
  }

  TString* tstr = ifile.Get<TString>("v_mask");
  if (tstr == nullptr) {
    throw std::runtime_error("Unable to read mask");
  }

  return std::bitset<3>(std::string(*tstr));
}

}  // namespace detail

ActsExamples::ScalingCalibrator::ScalingCalibrator(
    const std::filesystem::path& path)
    : m_calib_maps{::detail::readMaps(path)},
      m_mask{::detail::readMask(path)} {}

void ActsExamples::ScalingCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& /*gctx*/,
    Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const {
  Acts::SourceLink usl = trackState.getUncalibratedSourceLink();
  const IndexSourceLink& sourceLink = usl.get<IndexSourceLink>();

  assert((sourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  Acts::GeometryIdentifier mgid;
  mgid.setVolume(sourceLink.geometryId().volume() *
                 static_cast<Acts::GeometryIdentifier::Value>(m_mask[2]));
  mgid.setLayer(sourceLink.geometryId().layer() *
                static_cast<Acts::GeometryIdentifier::Value>(m_mask[1]));
  mgid.setSensitive(sourceLink.geometryId().sensitive() *
                    static_cast<Acts::GeometryIdentifier::Value>(m_mask[0]));
  const Cluster& cl = clusters->at(sourceLink.index());
  ConstantTuple ct = m_calib_maps.at(mgid).at(cl.sizeLoc0, cl.sizeLoc1);

  std::visit(
      [&](const auto& meas) {
        auto E = meas.expander();
        auto P = meas.projector();

        Acts::ActsVector<Acts::eBoundSize> fpar = E * meas.parameters();

        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * meas.covariance() * E.transpose();

        fpar[Acts::eBoundLoc0] += ct.x_offset;
        fpar[Acts::eBoundLoc1] += ct.y_offset;
        fcov(Acts::eBoundLoc0, Acts::eBoundLoc0) *= ct.x_scale;
        fcov(Acts::eBoundLoc1, Acts::eBoundLoc1) *= ct.y_scale;

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
