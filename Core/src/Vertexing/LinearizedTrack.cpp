// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/LinearizedTrack.hpp"

Acts::LinearizedTrack::LinearizedTrack()
{
	m_ParamsAtPCA.setZero();
	m_ParCovarianceAtPCA.setZero();
	m_LinPoint.setZero();
	m_PositionJacobian.setZero();
	m_MomentumJacobian.setZero();
	m_PositionAtPCA.setZero();
	m_MomentumAtPCA.setZero();
	m_ConstTerm.setZero();
}

Acts::LinearizedTrack::LinearizedTrack(const Acts::LinearizedTrack& other) :
m_ParamsAtPCA(other.m_ParamsAtPCA),
m_ParCovarianceAtPCA(other.m_ParCovarianceAtPCA),
m_LinPoint(other.m_LinPoint),
m_PositionJacobian(other.m_PositionJacobian),
m_MomentumJacobian(other.m_MomentumJacobian),
m_PositionAtPCA(other.m_PositionAtPCA),
m_MomentumAtPCA(other.m_MomentumAtPCA),
m_ConstTerm(other.m_ConstTerm) {}
 
Acts::LinearizedTrack::LinearizedTrack(const Acts::ActsVectorD<5>& paramsAtPCA,
								const Acts::ActsSymMatrixD<5>& parCovarianceAtPCA,
								const Acts::Vector3D& linPoint,
								const Acts::ActsMatrixD<5,3>& positionJacobian,
								const Acts::ActsMatrixD<5,3>& momentumJacobian,
								const Acts::Vector3D& positionAtPCA,
								const Acts::Vector3D& momentumAtPCA,
								const Acts::ActsVectorD<5>& constTerm) :
m_ParamsAtPCA(paramsAtPCA),
m_ParCovarianceAtPCA(parCovarianceAtPCA),
m_LinPoint(linPoint),
m_PositionJacobian(positionJacobian),
m_MomentumJacobian(momentumJacobian),
m_PositionAtPCA(positionAtPCA),
m_MomentumAtPCA(momentumAtPCA),
m_ConstTerm(constTerm) {}

Acts::LinearizedTrack& Acts::LinearizedTrack::operator= (const Acts::LinearizedTrack& other)
{
	if(this!=&other)
	{
		m_ParamsAtPCA 			= 		other.m_ParamsAtPCA;
		m_ParCovarianceAtPCA 	= 		other.m_ParCovarianceAtPCA;
		m_LinPoint 				= 		other.m_LinPoint;
		m_PositionJacobian		= 		other.m_PositionJacobian;
		m_MomentumJacobian		= 		other.m_MomentumJacobian;
		m_PositionAtPCA 		= 		other.m_PositionAtPCA;
		m_MomentumAtPCA 		= 		other.m_MomentumAtPCA;
		m_ConstTerm 			= 		other.m_ConstTerm;
	}
	return *this;
}


Acts::LinearizedTrack::~LinearizedTrack() {}