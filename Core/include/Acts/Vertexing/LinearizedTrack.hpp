// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once 

#include "Acts/Utilities/Definitions.hpp"

namespace Acts{

/**
 * @class LinearizedTrack.h
 *
 * Class for linear expansion of track parameters in vicinity of vertex
 *
 * The measurement equation is linearized in the following way: 
 *
 * F_k= D_k (x_k - x_0k) + E_k (p_k - p_0k) + F^0_k
 *
 * where F_k are the parameters at perigee nearest to the lin point, 
 * x_k is the position of the vertex, p_k the track momentum at the vertex,
 * and F^0_k is the constant term of expansion. D_k and E_k are matrices
 * of derivatives, denoted hereafter as "positionJacobian" and "momentumJacobian"
 * respectively.
 */

class LinearizedTrack
{
public:

	/// Default constructor
	LinearizedTrack();

	/// Copy constructor
	LinearizedTrack(const LinearizedTrack& other);

	/**
	 * Constructor taking perigee parameters and covariance matrix
	 * of track propagated to closest approach (PCA) of linearization point, 
	 * position and momentum Jacobian and const term.
	 */
	LinearizedTrack(const Acts::ActsVectorD<5>& paramsAtPCA,
					const Acts::ActsSymMatrixD<5>& parCovarianceAtPCA,
					const Acts::Vector3D& linPoint,
					const Acts::ActsMatrixD<5,3>& positionJacobian,
					const Acts::ActsMatrixD<5,3>& momentumJacobian,
					const Acts::Vector3D& positionAtPCA,
					const Acts::Vector3D& momentumAtPCA,
					const Acts::ActsVectorD<5>& constTerm);

	/// Assignment operator
	LinearizedTrack& operator= (const LinearizedTrack&);


	/// Default constructor
	virtual ~LinearizedTrack();

	const Acts::ActsVectorD<5>& parametersAtPCA() const
	{
		return m_ParamsAtPCA;
	}

	const Acts::ActsSymMatrixD<5>& covarianceAtPCA() const
	{
		return m_ParCovarianceAtPCA;
	}

	const Acts::Vector3D& linearizationPoint() const
	{
		return m_LinPoint;
	}

	const Acts::ActsMatrixD<5,3>& positionJacobian() const
	{
		return m_PositionJacobian;
	}

	const Acts::ActsMatrixD<5,3>& momentumJacobian() const
	{
		return m_MomentumJacobian;
	}

	const Acts::Vector3D& positionAtPCA() const
	{
		return m_PositionAtPCA;
	}

	const Acts::Vector3D& momentumAtPCA() const
	{
		return m_MomentumAtPCA;
	}

	const Acts::ActsVectorD<5>& constantTerm() const
	{
		return m_ConstTerm;
	}


private:

	 Acts::ActsVectorD<5> 		m_ParamsAtPCA;
	 Acts::ActsSymMatrixD<5> 	m_ParCovarianceAtPCA;
	 Acts::Vector3D 			m_LinPoint;
	 Acts::ActsMatrixD<5,3>		m_PositionJacobian;
	 Acts::ActsMatrixD<5,3>		m_MomentumJacobian;
	 Acts::Vector3D				m_PositionAtPCA;
	 Acts::Vector3D				m_MomentumAtPCA;
	 Acts::ActsVectorD<5>		m_ConstTerm;
};

} // namespace Acts
