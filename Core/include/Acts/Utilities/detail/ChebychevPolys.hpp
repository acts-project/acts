// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


/** @brief Acts implementation of the Chebychev polynomials & their first & second derivatives. 
 *         For a detailed explanation please refer to  the Wikipedia article 
 *              https://en.wikipedia.org/wiki/Chebyshev_polynomials
 * */
namespace Acts{

    namespace detail{
        /** @brief Returns the n-th Chebychev polynomial of first kind, T_{n}(x), evaluated at x. 
         *  @param k: Order of the polynomial
         *  @param x: Place of evaluation */
        constexpr double chebychevPolyTn(const unsigned k, const double x){
            switch (k) {
                case 0:
                    return 1;
                case 1:
                    return x;
                default:
                    return 2.*x*chebychevPolyTn(k-1, x) - chebychevPolyTn(k -2, x);
            }
        }
        /**  @brief Returns the n-th Chebychev polynomial of second kind, U_{n}(x), evaluated at x
          *  @param k: Order of the polynomial
          *  @param x: Place of evaluation */
        constexpr double chebychevPolyUn(const unsigned int k, const double x)  {
            switch (k) {
                case 0:
                    return 1.;
                    break;
                case 1:
                    return 2.*x;
                    break;
                default:
                    return 2.*chebychevPolyTn(k, x) + chebychevPolyUn(k-2, x);
            }
        }
    }
}