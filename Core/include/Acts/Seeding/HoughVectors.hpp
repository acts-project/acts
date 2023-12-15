// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <sstream>
#include <stdexcept>
#include <vector>

template <typename T>
class MultiIndexedVector2D {
 private:
  std::size_t m_d1, m_d2;
  std::vector<T> m_data;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

 public:
  MultiIndexedVector2D() : m_d1(0), m_d2(0) {
    m_logger = Acts::getDefaultLogger("HoughVectors", Acts::Logging::ERROR);
  }

  MultiIndexedVector2D(std::size_t d1, std::size_t d2, T const& t = T())
      : m_d1(d1), m_d2(d2), m_data(d1 * d2, t) {
    m_logger = Acts::getDefaultLogger("HoughVectors", Acts::Logging::ERROR);
  }

  std::size_t size(int dim) const {
    if (dim == 0) {
      return m_d1;
    }
    if (dim == 1) {
      return m_d2;
    } else {
      ACTS_ERROR("MultiIndexedVector2D: Argument to size() must be 0 or 1");
      return 0;
    }
  }

  void resize(std::size_t x1, std::size_t x2, T const& t = T()) {
    m_d1 = x1;
    m_d2 = x2;
    m_data.resize(x1 * x2, t);
  }

  T& operator()(std::size_t i, std::size_t j) {
    if (i >= m_d1 || j >= m_d2) {
      std::stringstream s;
      s << "MultiIndexedVector2D out of bounds: request (" << i << "," << j
        << ") size (" << m_d1 << "," << m_d2 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data[i * m_d2 + j];
  }
  inline T& operator()(std::pair<std::size_t, std::size_t> indices) {
    return operator()(indices.first, indices.second);
  }

  T const& operator()(std::size_t i, std::size_t j) const {
    if (i >= m_d1 || j >= m_d2) {
      std::stringstream s;
      s << "MultiIndexedVector2D out of bounds: request (" << i << "," << j
        << ") size (" << m_d1 << "," << m_d2 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data[i * m_d2 + j];
  }
  inline T const& operator()(
      std::pair<std::size_t, std::size_t> indices) const {
    return operator()(indices.first, indices.second);
  }

  T* operator[](std::size_t i) {
    if (i >= m_d1) {
      std::stringstream s;
      s << "MultiIndexedVector2D out of bounds: request " << i << " size ("
        << m_d1 << "," << m_d2 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data.data() + (i * m_d2);
  }

  const T* operator[](std::size_t i) const {
    if (i >= m_d1) {
      std::stringstream s;
      s << "MultiIndexedVector2D out of bounds: request " << i << " size ("
        << m_d1 << "," << m_d2 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data.data() + (i * m_d2);
  }

  void reset(const T& t) { m_data.assign(m_data.size(), t); }

  T* data() { return m_data.data(); }

  const T* data() const { return m_data.data(); }
};

template <typename T>
class MultiIndexedVector3D {
 private:
  std::size_t m_d1, m_d2, m_d3;
  std::vector<T> m_data;
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

 public:
  MultiIndexedVector3D() : m_d1(0), m_d2(0), m_d3(0) {
    m_logger = Acts::getDefaultLogger("HoughVectors", Acts::Logging::ERROR);
  }

  MultiIndexedVector3D(std::size_t d1, std::size_t d2, std::size_t d3,
                       T const& t = T())
      : m_d1(d1), m_d2(d2), m_d3(d3), m_data(d1 * d2 * d3, t) {
    m_logger = Acts::getDefaultLogger("HoughVectors", Acts::Logging::ERROR);
  }

  T& operator()(std::size_t i, std::size_t j, std::size_t k) {
    if (i >= m_d1 || j >= m_d2 || k >= m_d3) {
      std::stringstream s;
      s << "MultiIndexedVector3D out of bounds: request (" << i << "," << j
        << "," << k << ") size (" << m_d1 << "," << m_d2 << "," << m_d3 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data[i * m_d2 * m_d3 + j * m_d3 + k];
  }

  T const& operator()(std::size_t i, std::size_t j, std::size_t k) const {
    if (i >= m_d1 || j >= m_d2 || k >= m_d3) {
      std::stringstream s;
      s << "MultiIndexedVector3D out of bounds: request (" << i << "," << j
        << "," << k << ") size (" << m_d1 << "," << m_d2 << "," << m_d3 << ")";
      ACTS_ERROR(s.str());
    }
    return m_data[i * m_d2 * m_d3 + j * m_d3 + k];
  }

  void resize(std::size_t x1, std::size_t x2, std::size_t x3,
              T const& t = T()) {
    m_d1 = x1;
    m_d2 = x2;
    m_d3 = x3;
    m_data.resize(x1 * x2 * x3, t);
  }

  T* data() { return m_data.data(); }
};
