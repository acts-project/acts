// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"
#include "Acts/Definitions/Units.hpp"
using namespace Acts;
namespace Acts::Test{
    class TestSpacePoint{
            public:
                TestSpacePoint(const Vector3& pos,
                               const Vector3& dir,
                               const double radius,
                            const double time,
                            ActsSquareMatrix<3>& cov):
                            m_pos{pos},
                            m_dir{dir}, m_radius{radius}, m_time{time}, m_cov{cov}{}

                const Vector3&  localPosition() const {return m_pos; }
                const Vector3& sensorDirection() const { return m_dir;}

                const Vector3& stripPlaneNormal() const { return m_dir;}

                double driftRadius() const { return m_radius;}
                double time() const { return m_time;}
                const ActsSquareMatrix<3>& covariance() const { return m_cov;}
        private:
            Vector3 m_pos{Vector3::Zero()};
            Vector3 m_dir{Vector3::Zero()};
            double m_radius{0.};
            double m_time{0.};
            ActsSquareMatrix<3> m_cov{ActsSquareMatrix<3>::Zero()};
    };

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

BOOST_AUTO_TEST_CASE(WireResidualTest) {

    using namespace Acts::detail;
    using namespace Acts::UnitLiterals;
    using Line_t = StrawLineFitAuxiliaries::Line_t;
    using Vector = Line_t::Vector;
  

/// Set the line to be 45 degrees
    Line_t line{};
{
    using Pars_t = Line_t::ParamVector;
    using ParIdx = Line_t::ParIndices;
    Pars_t linePars{};
    linePars[ParIdx::phi] = 90 * 1_degree;
    linePars[ParIdx::theta] = 45 * 1_degree;
    line.updateParameters(linePars);
}

    const Vector wireDir = Vector::UnitZ();



}


}
