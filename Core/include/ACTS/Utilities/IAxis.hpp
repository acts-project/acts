// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts {

namespace detail {
  enum class AxisBoundaryType;
}

class IAxis
{

public:
  virtual bool
  isEquidistant() const = 0;
  virtual bool
  isVariable() const = 0;
  virtual detail::AxisBoundaryType
  getBoundaryType() const = 0;
  virtual std::vector<double>
  getBinEdges() const = 0;
  virtual double
  getMin() const = 0;
  virtual double
  getMax() const = 0;
  virtual size_t
  getNBins() const = 0;
};
}
