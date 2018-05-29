// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/variant/recursive_wrapper.hpp>
#include <boost/variant/variant_fwd.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Acts {

class variant_map;
class variant_vector;

/// Forward declaration for the variant data structure
/// using boost's own forward variant.
using variant_data = boost::variant<int,
                                    double,
                                    std::string,
                                    bool,
                                    boost::recursive_wrapper<variant_map>,
                                    boost::recursive_wrapper<variant_vector>>;
}