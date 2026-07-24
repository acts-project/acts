// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MultiAxisFactory.hpp"

#include <algorithm>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace Acts {

MultiAxisFactory::MultiAxisFactory(std::vector<AxisFactory> axisFactories)
    : m_axisFactories(std::move(axisFactories)) {
  if (m_axisFactories.empty()) {
    throw std::invalid_argument(
        "MultiAxisFactory: at least one axis description is required");
  }
}

std::size_t MultiAxisFactory::size() const {
  return m_axisFactories.size();
}

const AxisFactory& MultiAxisFactory::axisFactory(std::size_t i) const {
  return m_axisFactories.at(i);
}

const std::vector<AxisFactory>& MultiAxisFactory::axisFactories() const {
  return m_axisFactories;
}

bool MultiAxisFactory::isDeferred() const {
  return std::ranges::any_of(
      m_axisFactories, [](const AxisFactory& af) { return af.isDeferred(); });
}

std::vector<std::unique_ptr<IAxis>> MultiAxisFactory::toAxes() const {
  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(size());
  for (const AxisFactory& af : m_axisFactories) {
    axes.push_back(af.toAxis());
  }
  return axes;
}

std::vector<std::unique_ptr<IAxis>> MultiAxisFactory::toAxes(
    std::span<const AxisResolution> resolutions,
    std::span<const AxisDirection> directions) const {
  if (resolutions.size() != size()) {
    throw std::invalid_argument(
        "MultiAxisFactory: one resolution per axis is required");
  }
  if (!directions.empty() && directions.size() != size()) {
    throw std::invalid_argument(
        "MultiAxisFactory: one direction per axis is required");
  }
  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(size());
  for (std::size_t i = 0; i < size(); ++i) {
    std::optional<AxisDirection> direction;
    if (!directions.empty()) {
      direction = directions[i];
    }
    axes.push_back(m_axisFactories[i].toAxis(resolutions[i], direction));
  }
  return axes;
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::toMultiAxis() const {
  return makeMultiAxis(toAxes());
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::toMultiAxis(
    std::span<const AxisResolution> resolutions,
    std::span<const AxisDirection> directions) const {
  return makeMultiAxis(toAxes(resolutions, directions));
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::makeMultiAxis(
    const std::vector<std::unique_ptr<IAxis>>& axes) {
  switch (axes.size()) {
    case 1:
      return IMultiAxis::create(*axes[0]);
    case 2:
      return IMultiAxis::create(*axes[0], *axes[1]);
    case 3:
      return IMultiAxis::create(*axes[0], *axes[1], *axes[2]);
    default:
      throw std::domain_error(
          "MultiAxisFactory: multi-axes support at most 3 dimensions");
  }
}

std::string MultiAxisFactory::toString() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

std::ostream& operator<<(std::ostream& os,
                         const MultiAxisFactory& multiAxisFactory) {
  os << "MultiAxisFactory: " << multiAxisFactory.size() << " axes [";
  for (std::size_t i = 0; i < multiAxisFactory.size(); ++i) {
    os << (i > 0 ? "; " : "") << multiAxisFactory.axisFactory(i);
  }
  os << "]";
  return os;
}

}  // namespace Acts
