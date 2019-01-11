// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts{

template <typename BField>
class FullVertexFitter
{
public:

	struct Config
	{
		/// Magnetic field
		BField bField;
		/// Starting point of vertex fit
		Acts::Vector3D startingPoint;
		/// Maximum number of interations in fitter
		int maxIterations = 5;


		Config(BField bIn, Acts::Vector3D startP, int maxIter) :
			bField(std::move(bIn)),
			startingPoint(startP),
			maxIterations(maxIter)	
		{}

		/// Constructor with default number of iterations
		Config(BField bIn, Acts::Vector3D startP) :
			bField(std::move(bIn)),
			startingPoint(startP)
		{}			

		/// Constructor with default starting point
		Config(BField bIn, int maxIter) :
			bField(std::move(bIn)),
			startingPoint(Acts::Vector3D(0,0,0)),
			maxIterations(maxIter)
		{}

		/// Constructor with default number of iterations and starting point
		Config(BField bIn) :
			bField(std::move(bIn)),
			startingPoint(Acts::Vector3D(0,0,0))
		{}
	};

	/// Constructor with explicit config
	FullVertexFitter(const Config& cfg) : 
		m_cfg(cfg)
	{}

	/// Fit method, fitting vertex for provided tracks
	/// @param paramVector Vector of tracks to fit vertex to
	/// @return Fitted vertex
	Acts::Vertex fit(const std::vector<Acts::BoundParameters>& paramVector) const;

private:
	/// Configuration object
	Config m_cfg;
};

} // namespace Acts

#include "FullVertexFitter.ipp"