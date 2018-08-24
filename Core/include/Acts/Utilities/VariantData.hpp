// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// This file is a "Hello, world!" in C++ language by GCC for wandbox.
#include <boost/variant.hpp>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

namespace Acts {

/// Forward declarations for the recursive definition to work.
class variant_map;
class variant_vector;

/// This data type allows construction of a tree-like
/// data structure that contains varying data, in particular
/// int, double, strings, bools, maps and vectors. The latter
/// two can also contain variant_data (this allows trees)
/// The definition here uses boost's recursive_wrapper, because
/// @c variant_map and @c variant_vector have to be defined later on.
using variant_data = boost::variant<int,
                                    double,
                                    std::string,
                                    bool,
                                    boost::recursive_wrapper<variant_map>,
                                    boost::recursive_wrapper<variant_vector>>;

/// @brief Class which wraps an @c std::map<std::string, variant_data>
///
/// This class allows tree-like recursive data, and provides some
/// convenience methods to access keys. The key type is always @c std::string
class variant_map
{
public:
  /// Map type used by the class
  using map_t = std::map<std::string, variant_data>;

  using iterator       = map_t::iterator;
  using const_iterator = map_t::const_iterator;

  /// Constructor which takes a @c map_t and wraps it.
  /// @param src The source map, empty by default
  variant_map(map_t src = {}) : m_map(std::move(src)) {}

  /// Method to get the size of the map. Is forwarded to the
  /// contained map
  /// @return The size
  size_t
  size() const
  {
    return m_map.size();
  }

  /// Subscript operator to elements in the map
  /// @param key The key to access
  /// @return Value stored here. Note that this is @c variant_data.
  variant_data& operator[](const std::string& key) { return m_map[key]; }

  /// Method to count a key. Is forwarded
  /// @param key The key to count elements for
  /// @return Number of elements found
  size_t
  count(const std::string& key) const
  {
    return m_map.count(key);
  }

  /// Convenience method to access a key and transform it to a specific
  /// type. This will throw an error if the key is not found.
  /// @tparam T The type you want
  /// @param key The key to look up in the map
  template <class T>
  T&
  get(const std::string& key)
  {
    if (!m_map.count(key)) {
      throw std::out_of_range("variant_map key " + key + " not found");
    }
    return boost::get<T>(m_map.at(key));
  }

  /// Convenience method to access a key and transform it to a specific
  /// type. This will throw an error if the key is not found.
  /// @tparam T The type you want
  /// @param key The key to look up in the map
  template <class T>
  const T&
  get(const std::string& key) const
  {
    if (!m_map.count(key)) {
      throw std::out_of_range("variant_map key " + key + " not found");
    }
    return boost::get<T>(m_map.at(key));
  }

  /// Direct access to a key in the map, is forwarded
  /// @param key The key to look up
  /// @return value as @c variant_data
  const variant_data&
  at(const std::string& key) const
  {
    return m_map.at(key);
  }

  /// Direct access to a key in the map, is forwarded
  /// @param key The key to look up
  /// @return value as @c variant_data
  variant_data&
  at(const std::string& key)
  {
    return m_map.at(key);
  }

  /// Returns an iterator for the contained map, pointing to begin
  /// @return Iterator pointing at begin
  iterator
  begin()
  {
    return m_map.begin();
  }

  /// Returns an iterator for the contained map, pointing to end
  /// @return Iterator pointing at end
  iterator
  end()
  {
    return m_map.end();
  }

  /// Returns an const iterator for the contained map, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  begin() const
  {
    return m_map.begin();
  }

  /// Returns an const iterator for the contained map, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  end() const
  {
    return m_map.end();
  }

  /// Returns a const iterator for the contained map, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  cbegin()
  {
    return m_map.cbegin();
  }

  /// Returns a const iterator for the contained map, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  cend()
  {
    return m_map.cend();
  }

  /// Returns an const iterator for the contained map, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  cbegin() const
  {
    return m_map.cbegin();
  }

  /// Returns an const iterator for the contained map, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  cend() const
  {
    return m_map.cend();
  }

private:
  map_t m_map;
};

/// Class that wraps a vector of @c variant_data items
/// Provides some convenience methods for access.
class variant_vector
{
public:
  /// Vector type used by the class
  using vector_t       = std::vector<variant_data>;
  using iterator       = vector_t::iterator;
  using const_iterator = vector_t::const_iterator;

  /// Constructor for the class accepting an input vector
  /// @param src The source vector to wrap
  variant_vector(vector_t src = {}) : m_vector(std::move(src)) {}

  /// Method to get the current size of the wrapped vector
  /// @return The size
  size_t
  size() const
  {
    return m_vector.size();
  }

  /// Operator to get subscript access to an element.
  /// @param idx The index to access
  /// @return The value at idx, as @c variant_data
  variant_data& operator[](size_t idx) { return m_vector[idx]; }

  /// Method to access value at an index
  /// @param idx The index to access
  /// @return The value at idx, as @c variant_data
  variant_data&
  at(size_t idx)
  {
    return m_vector.at(idx);
  }

  /// Method to access value at an index
  /// @param idx The index to access
  /// @return The value at @c idx, as @c variant_data
  const variant_data&
  at(size_t idx) const
  {
    return m_vector.at(idx);
  }

  /// Method to append an element to the vector
  /// @param data @c variant_data value to add
  void
  push_back(const variant_data& data)
  {
    m_vector.push_back(data);
  }

  /// Convenience method which accesses a value and converts it
  /// to a provided type
  /// @tparam T the type you want
  /// @param idx The index to access
  /// @return The value at @c idx as @c T
  template <class T>
  T&
  get(const size_t& idx)
  {
    return boost::get<T>(m_vector.at(idx));
  }

  /// Convenience method which accesses a value and converts it
  /// to a provided type
  /// @tparam T the type you want
  /// @param idx The index to access
  /// @return The value at @c idx as @c T
  template <class T>
  const T&
  get(const size_t& idx) const
  {
    return boost::get<T>(m_vector.at(idx));
  }

  /// Returns an iterator for the contained vector, pointing to begin
  /// @return Iterator pointing at begin
  iterator
  begin()
  {
    return m_vector.begin();
  }

  /// Returns an iterator for the contained vector, pointing to end
  /// @return Iterator pointing at end
  iterator
  end()
  {
    return m_vector.end();
  }

  /// Returns an const iterator for the contained vector, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  begin() const
  {
    return m_vector.begin();
  }

  /// Returns an const iterator for the contained vector, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  end() const
  {
    return m_vector.end();
  }

  /// Returns an const iterator for the contained vector, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  cbegin()
  {
    return m_vector.cbegin();
  }

  /// Returns an const iterator for the contained vector, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  cend()
  {
    return m_vector.cend();
  }

  /// Returns an const iterator for the contained vector, pointing to begin
  /// @return Const iterator pointing at begin
  const_iterator
  cbegin() const
  {
    return m_vector.cbegin();
  }

  /// Returns an const iterator for the contained vector, pointing to end
  /// @return Const iterator pointing at end
  const_iterator
  cend() const
  {
    return m_vector.cend();
  }

private:
  vector_t m_vector;
};

/// This class implements a @c boost::static_vistor
/// which can turn a @c variant_data tree into a json string
class variant_json_visitor : public boost::static_visitor<>
{
public:
  /// Constructor for the visitor.
  /// @param pretty Whether the output is indented.
  variant_json_visitor(bool pretty = false)
    : m_json_str(std::stringstream()), m_pretty(pretty)
  {
    m_json_str << std::fixed << std::setprecision(30);
  }

  /// Visitor operator overload
  /// @param b The bool
  void
  operator()(const bool& b)
  {
    if (b) {
      m_json_str << "true";
    } else {
      m_json_str << "false";
    }
  }

  /// Visitor operator overload
  /// @param i The int
  void
  operator()(const int& i)
  {
    m_json_str << i;
  }

  /// Visitor operator overload
  /// @param d The double
  void
  operator()(const double& d)
  {
    m_json_str << d;
  }

  /// Visitor operator overload
  /// @param str The string
  void
  operator()(const std::string& str)
  {
    m_json_str << "\"" << str << "\"";
  }

  /// Visitor operator overload
  /// @param map The map
  void
  operator()(const variant_map& map)
  {
    m_json_str << "{";
    if (m_pretty) {
      m_json_str << std::endl;
    }

    size_t i = 0;
    m_depth += 1;
    for (const auto& entry : map) {
      indent();
      m_json_str << "\"" << entry.first << "\": ";
      boost::apply_visitor(*this, entry.second);

      if (i < map.size() - 1) {
        m_json_str << ", ";
        if (m_pretty) {
          m_json_str << std::endl;
        }
      }
      ++i;
    }
    m_depth -= 1;

    if (m_pretty) {
      m_json_str << std::endl;
    }
    indent();
    m_json_str << "}";
  }

  /// Visitor operator overload
  /// @param vec The vector
  void
  operator()(const variant_vector& vec)
  {
    m_json_str << "[";
    if (m_pretty) {
      m_json_str << std::endl;
    }
    size_t i = 0;
    m_depth += 1;
    for (const auto& entry : vec) {
      indent();
      boost::apply_visitor(*this, entry);
      if (i < vec.size() - 1) {
        m_json_str << ", ";
        if (m_pretty) {
          m_json_str << std::endl;
        }
      }
      ++i;
    }

    m_depth -= 1;
    if (m_pretty) {
      m_json_str << std::endl;
    }
    indent();
    m_json_str << "]";
  }

  /// Method to get the resulting json string
  /// @return The json string
  std::string
  str() const
  {
    return m_json_str.str();
  }

private:
  std::stringstream m_json_str;
  bool              m_pretty;
  size_t            m_depth = 0;

  void
  indent()
  {
    if (m_pretty) {
      for (size_t d = 0; d < m_depth; d++) {
        m_json_str << "  ";
      }
    }
  }
};

/// Function to apply the @c variant_json_visitor to a @c variant_data
/// instance. This turns the input into a json string and returns it.
/// @param data The input data to serialize
/// @param pretty Whether to indent the output
/// @return JSON string representation of the variant_data
inline std::string
to_json(const variant_data& data, bool pretty = false)
{
  variant_json_visitor jv(pretty);
  boost::apply_visitor(jv, data);
  return jv.str();
}

/// Operator overload for the stream operator. This allows
/// printing of @c variant_data instances.
/// @param os The outstream
/// @param data The variant_data
/// @return The stream given as @p os
inline std::ostream&
operator<<(std::ostream& os, const variant_data& data)
{
  os << to_json(data, true) << std::endl;
  return os;
}

/// Function to turn a @c Vector2D into @c variant_data
/// @param vec The vector to transform
/// @return The @c variant_data instance
inline variant_data
to_variant(const Vector2D& vec)
{
  using namespace std::string_literals;
  variant_map data;
  data["type"]    = "Vector2D"s;
  data["payload"] = variant_vector({vec[0], vec[1]});
  return data;
}

/// Function to turn a @c Transform3D into @c variant_data
/// @param trf The transform to convert
/// @return The @c variant_data instance
inline variant_data
to_variant(const Transform3D& trf)
{
  using namespace std::string_literals;
  variant_map data;
  data["type"] = "Transform3D"s;

  variant_map    payload;
  variant_vector matrix_data;
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < 4; j++) {
      matrix_data.push_back(trf(i, j));
    }
  }
  payload["data"] = matrix_data;

  data["payload"] = payload;
  return data;
}

/// Function to turn a @c ActsMatrix into @c variant_data
/// @param vec The vector to transform
/// @return The @c variant_data instance
inline variant_data
to_variant(const ActsMatrixD<4, 4>& matrix)
{
  using namespace std::string_literals;
  variant_map data;
  data["type"] = "Matrix4x4"s;

  variant_map payload;
  payload["cols"] = 4;
  payload["rows"] = 4;

  variant_vector matrix_data;
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < 4; j++) {
      matrix_data.push_back(matrix(i, j));
    }
  }
  payload["data"] = matrix_data;

  data["payload"] = payload;

  return data;
}

/// Function to produce an instance of type @c T from @c variant_data.
/// This is the unimplemented base template that is specialized
/// for various types.
/// @tparam T The type you want
/// @param vardata The data
/// @return The converted data as type @c T
template <typename T>
inline T
from_variant(const variant_data& vardata);

/// Function to produce a @c Transform3D from @c variant_data.
/// @param vardata The input @c variant_data
/// @return The converted @c Transform3D
template <>
inline Transform3D
from_variant<Transform3D>(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Must be variant_map");
  const variant_map& data = boost::get<variant_map>(vardata);
  throw_assert(data.get<std::string>("type") == "Transform3D",
               "Must be type Transform3D");

  const variant_map& payload = data.get<variant_map>("payload");

  const variant_vector& matrix_data = payload.get<variant_vector>("data");
  Transform3D           trf;
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < 4; j++) {

      size_t k     = i * 4 + j;
      double value = matrix_data.get<double>(k);
      trf(i, j) = value;
    }
  }

  return trf;
}

/// Function to produce a @c Vector2D from @c variant_data.
/// @param vardata The input @c variant_data
/// @return The converted @c Vector2D
template <>
inline Vector2D
from_variant<Vector2D>(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Must be variant_map");
  const variant_map& data = boost::get<variant_map>(vardata);
  throw_assert(data.get<std::string>("type") == "Vector2D",
               "Must be type Vector2D");

  const variant_vector& vector_data = data.get<variant_vector>("payload");

  Vector2D vec{vector_data.get<double>(0), vector_data.get<double>(1)};

  return vec;
}

}  // namespace Acts
