/*
 * eigen_migration
 * GeoPrimitivesCompare.h
 *
 *  Created on: Mar 25, 2014
 *      Author: rbianchi <Riccardo.Maria.Bianchi@cern.ch>
 *
 */

#ifndef GEOPRIMITIVESCOMPARE_H_
#define GEOPRIMITIVESCOMPARE_H_

#include "GeoPrimitives/GeoPrimitives.h"

namespace Amg {


/*
 * comparison of two Vector3D,
 * needed for the definition of a std::set<Amg::Vector3D>
 *
 * originals for CLHEP::HepVector3D was in CLHEP/Vector/src/SpaceVector.cc:
 *   int Hep3Vector::compare (const Hep3Vector & v) const;
 *   bool Hep3Vector::operator > (const Hep3Vector & v) const;
 *   bool Hep3Vector::operator < (const Hep3Vector & v) const;
 *
 */

struct Vector3DComparer
{
	bool operator()(Amg::Vector3D const &v1, Amg::Vector3D const &v2)
	{
		if ( v1.z() > v2.z() ) {
			return false;
		} else if ( v1.z() < v2.z() ) {
			return true;
		} else if ( v1.y() > v2.y() ) {
			return false;
		} else if ( v1.y() < v2.y() ) {
			return true;
		} else if ( v1.x() > v2.x() ) {
			return false;
		} else if ( v1.x() < v2.x() ) {
			return true;
		} else {
			return false;
		}
	}
};


struct VectorVector3DComparer
{
	bool operator()(std::vector<Amg::Vector3D> const &vv1, std::vector<Amg::Vector3D> const &vv2)
	{
		std::vector<Amg::Vector3D>::const_iterator v1;
		std::vector<Amg::Vector3D>::const_iterator v2;

		for (v1 = vv1.begin(); v1 != vv1.end(); ++v1) {
			for (v2 = vv2.begin(); v2 != vv2.end(); ++v2) {
				if ( (*v1).z() > (*v2).z() ) {
					return false;
				} else if ( (*v1).z() < (*v2).z() ) {
					return true;
				} else if ( (*v1).y() > (*v2).y() ) {
					return false;
				} else if ( (*v1).y() < (*v2).y() ) {
					return true;
				} else if ( (*v1).x() > (*v2).x() ) {
					return false;
				} else if ( (*v1).x() < (*v2).x() ) {
					return true;
				} else {
					return false;
				}
			}
		}
		return false; // you should never get here
	}
};





} // end namespace

#endif /* GEOPRIMITIVESCOMPARE_H_ */
