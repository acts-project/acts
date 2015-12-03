/*
 * eigen_migration
 * EulerAnglesHelpers.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rbianchi <Riccardo.Maria.Bianchi@cern.ch>
 *
 */

#ifndef _GEOPRIMITIVES_EULERANGLESHELPERS_H
#define _GEOPRIMITIVES_EULERANGLESHELPERS_H

#include "GeoPrimitives/GeoPrimitives.h"

#include "CLHEPtoEigenEulerAnglesConverters.h"

namespace Amg {

/**
 * Get the equivalents to
 * CLHEP Phi, Theta, Psi Euler angles
 *
 * phi = vector[0]
 * theta = vector[1]
 * psi = vector[2]
 *
 * N.B.
 * if "convention = 0" --> "Z-X-Z" convention ==> DEFAULT!!
 * if "convention = 1" --> "Z-Y-Z" convention
 *
 * N.B.!!
 * for normal usage, use the default notation (simply leave it empty, or use convention=0),
 * or, alternatively, be sure to use the same convention in both setPhi() and getPhiThetaPsi().
 *
 *
 */
inline Amg::Vector3D getPhiThetaPsi(Amg::RotationMatrix3D mat, int convention = 0)
{
	double phi;
	double theta;
	double psi;

	Amg::Vector3D ea; // euler angles vector

	/**
	 * we extract the Euler Angles
	 * from the Eigen matrix,
	 *
	 */
	// using Z-X-Z convention: (2,0,2) == "Z,X,Z"
	if (convention == 0) {
		ea = mat.eulerAngles(2, 0, 2); // (using Z-X-Z convention) // DEFAULT
	}
	// using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
	else {
		ea = mat.eulerAngles(2, 1, 2); // (using Z-Y-Z convention)
	}

	// convert the values from Eigen convention to CLHEP convention
	ea = convert_EigenEulerAngles_to_CLHEPPhiThetaPsi( ea, convention );

	phi = ea[0];
	theta = ea[1];
	psi = ea[2];

	Amg::Vector3D phiThetaPsi_angles;
	phiThetaPsi_angles(0) = phi;
	phiThetaPsi_angles(1) = theta;
	phiThetaPsi_angles(2) = psi;

	return phiThetaPsi_angles;

} // end of getPhiThetaPsi()







/* Set the Phi angle of a rotation matrix,
 * leaving Theta and Psi unaltered.
 *
 * Actually it returns a new rotation matrix
 * built with the new Phi angle, and with the
 * Theta and Psi angles taken from the original matrix.
 *
 * N.B.
 * if "convention = 0" --> "Z-X-Z" convention ==> DEFAULT!!
 * if "convention = 1" --> "Z-Y-Z" convention
 *
 * N.B.!!
 * for normal usage, use the default notation (simply leave it empty, or use convention=0),
 * or, alternatively, be sure to use the same convention in both setPhi() and getPhiThetaPsi().
 *
 *
 */
inline Amg::RotationMatrix3D setPhi(Amg::RotationMatrix3D mat, double angle, int convention = 0)
{


	Amg::Vector3D phi_theta_psi;


	// using Z-X-Z convention: (2,0,2) == "Z,X,Z". ==>  DEFAULT!
	if (convention == 0) {
		phi_theta_psi = getPhiThetaPsi(mat, 0); // using Z-X-Z ((2,0,2) convention // DEFAULT
	}
	else { // using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
		phi_theta_psi = getPhiThetaPsi(mat, 1); // using Z-Y-Z ((2,1,2) convention
	}

	double phi = phi_theta_psi(0);
	double theta = phi_theta_psi(1);
	double psi = phi_theta_psi(2);


	/*
	 * set a new Phi angle, from user's settings
	 */
	phi = angle;

	/*
	 * get Eigen Euler angles from CLEHP style Phi, Theta, Psi
	 */
	Amg::Vector3D clhep_phiThetaPsi(phi, theta, psi); // a vector with CLHEP angles
	Amg::Vector3D eigen_euler_angles;
	// converting the CLHEP angles to Eigen angles
	if (convention == 0) { // using Z-X-Z convention: (2,0,2) == "Z,X,Z". ==>  DEFAULT!
		eigen_euler_angles = convert_CLHEPPhiThetaPsi_to_EigenEulerAngles(clhep_phiThetaPsi, 0); // using Z-X-Z ((2,0,2) convention // DEFAULT
	}
	else { // using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
		eigen_euler_angles = convert_CLHEPPhiThetaPsi_to_EigenEulerAngles(clhep_phiThetaPsi, 1); // using Z-Y-Z ((2,1,2) convention
	}

	/*
	 * create a new rotation matrix from AngleAxis
	 *
	 * N.B.!!
	 * to match CLHEP results,
	 * we have to invert the order of the rotation operations,
	 * compared to the order in CLHEP.
	 * The matrix here below is equal to:
	 * ---
	 * CLHEP::HepRotation localRot;
	 * localRot.rotateZ(angC);
	 * localRot.rotateY(angB);
	 * localRot.rotateZ(angA);
	 * ---
	 *
	 */

	if (convention == 0) { // using Z-X-Z convention: (2,0,2) == "Z,X,Z". ==>  DEFAULT!
		mat = Amg::AngleAxis3D(eigen_euler_angles(0), Amg::Vector3D::UnitZ())
		* Amg::AngleAxis3D(eigen_euler_angles(1), Amg::Vector3D::UnitX())
		* Amg::AngleAxis3D(eigen_euler_angles(2), Amg::Vector3D::UnitZ());
	}
	else { // using Z-Y-Z convention: (2,1,2) == "Z,Y,Z"
		mat = Amg::AngleAxis3D(eigen_euler_angles(0), Amg::Vector3D::UnitZ())
		* Amg::AngleAxis3D(eigen_euler_angles(1), Amg::Vector3D::UnitY())
		* Amg::AngleAxis3D(eigen_euler_angles(2), Amg::Vector3D::UnitZ());
	}

	return mat;

} // end of SetPhi()




} // end of namespace Amg


#endif

