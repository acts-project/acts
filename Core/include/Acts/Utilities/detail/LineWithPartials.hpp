// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"


namespace Acts::detail {
    template <std::floating_point T>
    class LineWithPartials {
        public:
        LineWithPartials() = default;

        enum ParIndices : unsigned{
            x0= 0,
            y0 =1,
            theta= 2,
            phi =3,
            nPars =4
        };
        using Vector = Eigen::Matrix<T, 3, 1>;
        using ParamVector = Eigen::Matrix<T, nPars, 1>;
        
        void updateParameters(const ParamVector& newPars);
        
        const Vector& position() const; 
        const Vector& direction() const;
        const Vector& gradient(const unsigned param) const;
        const Vector& hessian(const unsigned param1, const unsigned param2) const;
        private:
        Vector m_pos{Vector::Zero()};
        Vector m_dir{Vector::Zero()};
        static constexpr std::size_t s_nPars = ParIndices::nPars;
        std::array<Vector, s_nPars> m_gradient{filledArray<Vector, s_nPars>(Vector::Zero())};
        
        std::array<Vector, sumUpToN(s_nPars)> m_hessian{filledArray<Vector, sumUpToN(s_nPars)>(Vector::Zero())};   
    };
}
#include "Acts/Utilities/detail/LineWithPartials.ipp"