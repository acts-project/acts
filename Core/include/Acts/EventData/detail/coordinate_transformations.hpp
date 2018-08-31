// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// Acts includes
#include <cmath>
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#ifdef ACTS_COORDINATE_TRANSFORM_PLUGIN

#include ACTS_COORDINATE_TRANSFORM_PLUGIN

#else

namespace Acts {
/// @cond detail
namespace detail {
  /**
   * @brief helper structure summarizing coordinate transformations
   */
  struct coordinate_transformation
  {
    using ParVector_t = ActsVector<ParValue_t, Acts::NGlobalPars>;

    static ActsVectorD<3>
    parameters2globalPosition(const ParVector_t& pars, const Surface& s)
    {
      ActsVectorD<3> globalPosition;
      s.localToGlobal(ActsVectorD<2>(pars(Acts::eLOC_0), pars(Acts::eLOC_1)),
                      parameters2globalMomentum(pars),
                      globalPosition);
      return globalPosition;
    }

    static ActsVectorD<3>
    parameters2globalMomentum(const ParVector_t& pars)
    {
      ActsVectorD<3> momentum;
      double         p     = std::abs(1. / pars(Acts::eQOP));
      double         phi   = pars(Acts::ePHI);
      double         theta = pars(Acts::eTHETA);
      momentum << p * sin(theta) * cos(phi), p * sin(theta) * sin(phi),
          p * cos(theta);

      return momentum;
    }

    static ParVector_t
    global2curvilinear(const ActsVectorD<3>& /*pos*/,
                       const ActsVectorD<3>& mom,
                       double                charge)
    {
      using VectorHelpers::phi;
      using VectorHelpers::theta;
      ParVector_t parameters;
      parameters << 0, 0, phi(mom), theta(mom),
          ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm();

      return parameters;
    }

    static ParVector_t
    global2parameters(const ActsVectorD<3>& pos,
                      const ActsVectorD<3>& mom,
                      double                charge,
                      const Surface&        s)
    {
      using VectorHelpers::phi;
      using VectorHelpers::theta;
      ActsVectorD<2> localPosition;
      s.globalToLocal(pos, mom, localPosition);
      ParVector_t result;
      result << localPosition(0), localPosition(1), phi(mom), theta(mom),
          ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm();
      return result;
    }

    static double
    parameters2charge(const ParVector_t& pars)
    {
      return (pars(Acts::eQOP) > 0) ? 1. : -1.;
    }
  };
}  // namespace detail
/// @endcond
}  // namespace Acts
#endif  // ACTS_COORDINATE_TRANSFORM_PLUGIN
