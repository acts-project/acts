// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace Acts {

namespace detail {

  template <bool equidistant>
  class Axis;

  typedef Axis<true>  EquidistantAxis;
  typedef Axis<false> VariableAxis;

  template <>
  class Axis<true>
  {
  public:
    Axis(double xmin, double xmax, size_t bins)
      : m_min(xmin), m_max(xmax), m_width((xmax - xmin) / bins), m_bins(bins)
    {
    }

    size_t
    getBin(double x) const
    {
      if (x <= getMin()) return 0u;
      if (x >= getMax()) return getNBins() - 1;

      return std::floor((x - getMin()) / getBinWidth());
    }

    double
    getBinWidth() const
    {
      return m_width;
    }

    double
    getMax() const
    {
      return m_max;
    }

    double
    getMin() const
    {
      return m_min;
    }

    size_t
    getNBins() const
    {
      return m_bins + 1;
    }

  private:
    double m_min;
    double m_max;
    double m_width;
    size_t m_bins;
  };

  template <>
  class Axis<false>
  {
  public:
    Axis(std::vector<double> binEdges) : m_binEdges(std::move(binEdges)) {}

    size_t
    getBin(double x) const
    {
      if (x <= getMin()) return 0u;
      if (x >= getMax()) return getNBins() - 1;

      const auto it
          = std::upper_bound(std::begin(m_binEdges), std::end(m_binEdges), x);
      return std::distance(std::begin(m_binEdges), it) - 1;
    }

    double
    getMax() const
    {
      return m_binEdges.back();
    }

    double
    getMin() const
    {
      return m_binEdges.front();
    }

    size_t
    getNBins() const
    {
      return m_binEdges.size();
    }

  private:
    std::vector<double> m_binEdges;
  };
}  // namespace detail

}  // namespace Acts
