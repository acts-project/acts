// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_UTILITIES_VARIANTDATA_H
#define ACTS_UTILITIES_VARIANTDATA_H 1

// This file is a "Hello, world!" in C++ language by GCC for wandbox.
#include <boost/variant.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ThrowAssert.hpp"

namespace Acts {

// typedef boost::make_recursive_variant<
// int,
// double,
// std::string,
// bool,
// std::map<std::string, boost::recursive_variant_>,
// std::vector<boost::recursive_variant_>>::type variant_data;

// typedef std::map<std::string, gen_structure> variant_map;
// typedef std::vector<gen_structure>           variant_vector;

class variant_map;
class variant_vector;

using variant_data = boost::variant<int,
                                    double,
                                    std::string,
                                    bool,
                                    boost::recursive_wrapper<variant_map>,
                                    boost::recursive_wrapper<variant_vector>>;

class variant_map
{
public:
  using map_t          = std::map<std::string, variant_data>;
  using iterator       = map_t::iterator;
  using const_iterator = map_t::const_iterator;

  variant_map(map_t src = {}) : m_map(src) {}

  size_t
  size() const
  {
    return m_map.size();
  }
  variant_data& operator[](std::string key) { return m_map[key]; }

  size_t
  count(const std::string& key) const
  {
    return m_map.count(key);
  }

  template <class T>
  T&
  get(const std::string& key)
  {
    return boost::get<T>(m_map.at(key));
  }

  template <class T>
  const T&
  get(const std::string& key) const
  {
    return boost::get<T>(m_map.at(key));
  }

  iterator
  begin()
  {
    return m_map.begin();
  }
  iterator
  end()
  {
    return m_map.end();
  }
  const_iterator
  begin() const
  {
    return m_map.begin();
  }
  const_iterator
  end() const
  {
    return m_map.end();
  }

private:
  map_t m_map;
};

class variant_vector
{
public:
  using vector_t       = std::vector<variant_data>;
  using iterator       = vector_t::iterator;
  using const_iterator = vector_t::const_iterator;

  variant_vector(vector_t src = {}) : m_vector(src) {}

  size_t
  size() const
  {
    return m_vector.size();
  }
  variant_data& operator[](size_t idx) { return m_vector[idx]; }

  variant_data&
  at(size_t idx)
  {
    return m_vector.at(idx);
  }

  const variant_data&
  at(size_t idx) const
  {
    return m_vector.at(idx);
  }

  void
  push_back(variant_data&& data)
  {
    m_vector.push_back(data);
  }
  
  template <class T>
  T&
  get(const size_t& idx)
  {
    return boost::get<T>(m_vector.at(idx));
  }

  template <class T>
  const T&
  get(const size_t& idx) const
  {
    return boost::get<T>(m_vector.at(idx));
  }

  iterator
  begin()
  {
    return m_vector.begin();
  }
  iterator
  end()
  {
    return m_vector.end();
  }
  const_iterator
  begin() const
  {
    return m_vector.begin();
  }
  const_iterator
  end() const
  {
    return m_vector.end();
  }

private:
  vector_t m_vector;
};

class variant_json_visitor : public boost::static_visitor<>
{
public:
  variant_json_visitor(bool pretty_ = false)
    : json_str(std::stringstream()), pretty(pretty_)
  {
  }

  void
  operator()(const bool& b)
  {
    if (b) {
      json_str << "true";
    } else {
      json_str << "false";
    }
  }

  void
  operator()(const int& i)
  {
    json_str << i;
  }

  void
  operator()(const double& d)
  {
    json_str << d;
  }

  void
  operator()(const std::string& str)
  {
    json_str << "\"" << str << "\"";
  }

  void
  operator()(const variant_map& map)
  {
    json_str << "{";
    if (pretty) json_str << std::endl;

    size_t i = 0;
    depth += 1;
    for (const auto& entry : map) {
      indent();
      json_str << "\"" << entry.first << "\": ";
      boost::apply_visitor(*this, entry.second);

      if (i < map.size() - 1) {
        json_str << ", ";
        if (pretty) json_str << std::endl;
      }
      ++i;
    }
    depth -= 1;

    if (pretty) json_str << std::endl;
    indent();
    json_str << "}";
  }

  void
  operator()(const variant_vector& vec)
  {
    json_str << "[";
    if (pretty) json_str << std::endl;
    size_t i = 0;
    depth += 1;
    for (const auto& entry : vec) {
      indent();
      boost::apply_visitor(*this, entry);
      if (i < vec.size() - 1) {
        json_str << ", ";
        if (pretty) json_str << std::endl;
      }
      ++i;
    }

    depth -= 1;
    if (pretty) json_str << std::endl;
    indent();
    json_str << "]";
  }

  std::string
  str() const
  {
    return json_str.str();
  }

private:
  std::stringstream json_str;
  bool              pretty;
  size_t            depth = 0;

  void
  indent()
  {
    if (pretty) {
      for (size_t d = 0; d < depth; d++) {
        json_str << "  ";
      }
    }
  }
};

inline std::string
to_json(const variant_data& data, bool pretty = false)
{
  variant_json_visitor jv(pretty);
  boost::apply_visitor(jv, data);
  return jv.str();
}

inline std::ostream&
operator<<(std::ostream& os, const variant_data& data)
{
  os << to_json(data, true) << std::endl;
  return os;
}

inline variant_map
to_variant(const Vector2D& vec)
{
  using namespace std::string_literals;
  variant_map data;
  data["type"]    = "Vector2D"s;
  data["payload"] = variant_vector({vec[0], vec[1]});
  return data;
}

inline variant_map
to_variant(const Transform3D& trf)
{
  using namespace std::string_literals;
  variant_map data;
  data["type"] = "Transform3D"s;

  variant_map payload;
  // payload["matrix"] = to_variant(trf.matrix());
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

inline variant_map
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

template<typename T>
inline T
from_variant(const variant_data& data_);

template<>
inline Transform3D
from_variant<Transform3D>(const variant_data& data_)
{
  throw_assert(data_.which() == 4, "Must be variant_map");
  const variant_map& data = boost::get<variant_map>(data_);
  throw_assert(data.get<std::string>("type") == "Transform3D",
               "Must be type Transform3D");

  const variant_map &payload = data.get<variant_map>("payload");

  const variant_vector &matrix_data = payload.get<variant_vector>("data");
  Transform3D trf;
  for(size_t i=0;i<4;i++) {
    for(size_t j=0;j<4;j++) {

      size_t k = i*4+j;
      double value = matrix_data.get<double>(k);
      trf(i, j) = value;
    }
  }

  return trf;
}

template<>
inline Vector2D
from_variant<Vector2D>(const variant_data& data_)
{
  throw_assert(data_.which() == 4, "Must be variant_map");
  const variant_map& data = boost::get<variant_map>(data_);
  throw_assert(data.get<std::string>("type") == "Vector2D",
               "Must be type Vector2D");

  const variant_vector &vector_data = data.get<variant_vector>("payload");

  Vector2D vec;
  for(size_t i=0;i<2;i++) {
    vec[i] = vector_data.get<double>(i);
  }

  return vec;
}

/*int
main()
{
  using namespace std::string_literals;

  auto data       = gen_map();
  data["kint"]    = 5;
  data["kdouble"] = 5.3;
  data["kstring"] = "whoop"s;
  data["kbool"]   = true;
  data["kvector"] = gen_vector({"a"s, 8, "c"s});
  data["kmap"]    = gen_map(
      {{"h", 17}, {"l", "hurz"s}, {"b", false}, {"w", gen_map({{"m", -1}})}});

  gen_structure root = data;

  std::cout << std::endl << std::endl;
  std::cout << root << std::endl;
}*/

}  // namespace Acts

#endif
