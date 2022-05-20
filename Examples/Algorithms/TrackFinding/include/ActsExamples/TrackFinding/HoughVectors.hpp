/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef HTTVECTORS_H
#define HTTVECTORS_H

/**
 * @file HTTVectors.h
 * @author Riley Xu - rixu@cern.ch
 * @date Jan 27 2020
 * @brief Defines several vector wrappers for homogenous multi-dimensional vectors, declared as 1D arrays for locality.
 */


#include <vector>
#include <stdexcept>

// #define DEBUG_HTTVECTORS

#ifdef DEBUG_HTTVECTORS
#include <sstream>
#endif


template <typename T>
class vector2D
{
 private:

  size_t m_d1,m_d2;
  std::vector<T> m_data;

 public:

 vector2D() :
  m_d1(0), m_d2(0)
    {}

 vector2D(size_t d1, size_t d2, T const & t=T()) :
  m_d1(d1), m_d2(d2), m_data(d1*d2, t)
    {}

  size_t size(int dim) const
  {
    if (dim == 0) return m_d1;
    if (dim == 1) return m_d2;
    else throw std::range_error("vector2D: Argument to size() must be 0 or 1");
  }

  void resize(size_t x1, size_t x2, T const & t=T())
  {
    m_d1 = x1;
    m_d2 = x2;
    m_data.resize(x1*x2, t);
  }

  T & operator()(size_t i, size_t j)
    {
#ifdef DEBUG_HTTVECTORS
      if (i >= m_d1 || j >= m_d2)
        {
	  std::stringstream s;
	  s << "vector2D out of bounds: request (" << i << "," << j << ") size (" << m_d1 << "," << m_d2 << ")";
	  throw std::range_error(s.str());
        }
#endif
      return m_data[i*m_d2 + j];
    }

  T const & operator()(size_t i, size_t j) const
  {
#ifdef DEBUG_HTTVECTORS
    if (i >= m_d1 || j >= m_d2)
      {
	std::stringstream s;
	s << "vector2D out of bounds: request (" << i << "," << j << ") size (" << m_d1 << "," << m_d2 << ")";
	throw std::range_error(s.str());
      }
#endif
    return m_data[i*m_d2 + j];
  }

  T* operator[](size_t i)
  {
#ifdef DEBUG_HTTVECTORS
    if (i >= m_d1)
      {
	std::stringstream s;
	s << "vector2D out of bounds: request " << i << " size (" << m_d1 << "," << m_d2 << ")";
	throw std::range_error(s.str());
      }
#endif
    return m_data.data() + (i * m_d2);
  }

  const T* operator[](size_t i) const
  {
#ifdef DEBUG_HTTVECTORS
    if (i >= m_d1)
      {
	std::stringstream s;
	s << "vector2D out of bounds: request " << i << " size (" << m_d1 << "," << m_d2 << ")";
	throw std::range_error(s.str());
      }
#endif
    return m_data.data() + (i * m_d2);
  }

  T* data()
  {
    return m_data.data();
  }

  const T* data() const
  {
    return m_data.data();
  }
};


template <typename T>
class vector3D
{
 private:

  size_t m_d1,m_d2,m_d3;
  std::vector<T> m_data;

 public:

 vector3D() :
  m_d1(0), m_d2(0), m_d3(0)
    {}

 vector3D(size_t d1, size_t d2, size_t d3, T const & t=T()) :
  m_d1(d1), m_d2(d2), m_d3(d3), m_data(d1*d2*d3, t)
    {}

  T & operator()(size_t i, size_t j, size_t k)
    {
#ifdef DEBUG_HTTVECTORS
      if (i >= m_d1 || j >= m_d2 || k >= m_d3)
        {
	  std::stringstream s;
	  s << "vector3D out of bounds: request (" << i << "," << j << "," << k << ") size (" << m_d1 << "," << m_d2 << "," << m_d3 << ")";
	  throw std::range_error(s.str());
        }
#endif
      return m_data[i*m_d2*m_d3 + j*m_d3 + k];
    }

  T const & operator()(size_t i, size_t j, size_t k) const
  {
#ifdef DEBUG_HTTVECTORS
    if (i >= m_d1 || j >= m_d2 || k >= m_d3)
      {
	std::stringstream s;
	s << "vector3D out of bounds: request (" << i << "," << j << "," << k << ") size (" << m_d1 << "," << m_d2 << "," << m_d3 << ")";
	throw std::range_error(s.str());
      }
#endif
    return m_data[i*m_d2*m_d3 + j*m_d3 + k];
  }

  void resize(size_t x1, size_t x2, size_t x3, T const & t=T())
  {
    m_d1 = x1;
    m_d2 = x2;
    m_d3 = x3;
    m_data.resize(x1*x2*x3, t);
  }

  T* data()
  {
    return m_data.data();
  }
};


#endif
