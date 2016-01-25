/*
 * eigen_migration
 * main.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: rbianchi <Riccardo.Maria.Bianchi@cern.ch>
 *
 */

#ifndef _GEOPRIMITIVES_CLHEPTOEIGENEULERANGLESCONVERTERS_H
#define _GEOPRIMITIVES_CLHEPTOEIGENEULERANGLESCONVERTERS_H

#include "GeoPrimitives/GeoPrimitives.h"

#include <math.h> // for M_PI definition



namespace Amg {


/**
 * Convert CLEHP Phi,Theta,Psi angles
 * to Eigen euler angles using Z-X-Z convention
 *
 * N.B.
 * if "convention = 0" --> "Z-X-Z" convention ==> DEFAULT!!
 * if "convention = 1" --> "Z-Y-Z"
 */
inline Amg::Vector3D convert_CLHEPPhiThetaPsi_to_EigenEulerAngles(Amg::Vector3D clhep_angles, int convention = 0)
{
	Amg::Vector3D eigen_angles;

	// using Z-X-Z convention (CLHEP): (2,0,2) == "Z,X,Z". // DEFAULT
	if (convention == 0) {
		eigen_angles(2) = -clhep_angles(0); // Phi-->Z
		eigen_angles(1) = -clhep_angles(1); // Theta-->X
		eigen_angles(0) = -clhep_angles(2); // Psi-->Z
	}
	// using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
	else {
		eigen_angles(0) = -0.5 * clhep_angles(0); // Phi-->Z
		eigen_angles(1) = clhep_angles(1); // Theta-->Y
		eigen_angles(2) = clhep_angles(2); // Psi-->Z
	}

	return eigen_angles;
}

/**
 * Convert Eigen euler angles
 * to CLEHP Phi,Theta,Psi angles
 *
 * N.B.
 * if "convention = 0" --> "Z-X-Z" convention ==> DEFAULT!!
 * if "convention = 1" --> "Z-Y-Z" convention
 */
inline Amg::Vector3D convert_EigenEulerAngles_to_CLHEPPhiThetaPsi(Amg::Vector3D eigen_angles, int convention = 0)
{
	Amg::Vector3D clhep_angles;

	double phi;
	double theta;
	double psi;

	/**
	 * Note:
	 * as explained in:  eigen / Eigen / src / Geometry / EulerAngles.h
	 * the returned angles are in the ranges [0:pi]x[-pi:pi]x[-pi:pi]
	 *
	 * (source here:
	 * https://bitbucket.org/eigen/eigen/src/42e011583bceb055a43fa688622e828fbbabf818/Eigen/src/Geometry/EulerAngles.h)
	 *
	 *
	 * N.B.!!
	 * CLHEP's Phi, Theta, Psi correspond to
	 * eulerAngles[2], [1] and [0] respectively,
	 * with the sign inverted.
	 */
	// using Z-X-Z convention: (2,0,2) == "Z,X,Z".  // DEFAULT
	if (convention == 0) {
		phi   = -eigen_angles(2); // Z-->Phi
		theta = -eigen_angles(1); // X-->Theta
		psi   = -eigen_angles(0); // Z-->Psi
	}
	// using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
	else {
		phi   = -2 * eigen_angles(0); // Z-->Phi
		theta = eigen_angles(1); // Y-->Theta
		psi   = eigen_angles(2); // Z-->Psi
	}

	clhep_angles(0) = phi;
	clhep_angles(1) = theta;
	clhep_angles(2) = psi;

	return clhep_angles;

}


} // end of namespace Amg


#endif





