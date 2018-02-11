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

  void
  push_back(variant_data&& data) {
    m_vector.push_back(data);
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

inline
std::string
to_json(const variant_data& data, bool pretty = false)
{
  variant_json_visitor jv(pretty);
  boost::apply_visitor(jv, data);
  return jv.str();
}

inline
std::ostream&
operator<<(std::ostream& os, const variant_data& data)
{
  os << to_json(data, true) << std::endl;
  return os;
}

inline
variant_map
to_variant(const Vector2D &vec) {
  using namespace std::string_literals;
  variant_map data;
  data["type"] = "Vector2D"s;
  data["payload"] = variant_vector({vec[0], vec[1]});
  return data;
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
