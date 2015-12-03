///////////////////////////////////////////////////////////////////
// GeoPrimitivesHelpers.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef GEOPRIMITIVES_GEOPRIMITIVESHELPERS_H
#define GEOPRIMITIVES_GEOPRIMITIVESHELPERS_H

#include "GeoPrimitives/GeoPrimitives.h"
#include "GeoPrimitives/GeoPrimitivesCompare.h"

#include "cmath"

#include <vector>
#include <set>
#include <iostream>


/** Geometry primitives helper functions
 @author  Niels van Eldik
 @author  Robert Johannes Langenberg <robert.langenberg@cern.ch>
 @author Andreas Salzburger <Andreas.Salzburger@cern.ch>

 */

namespace Amg {



typedef std::set<Amg::Vector3D, Vector3DComparer> SetVector3D;
typedef std::set< std::vector< Amg::Vector3D>, VectorVector3DComparer> SetVectorVector3D;



/** calculates the opening angle between two vectors */
inline double angle(const Amg::Vector3D& v1, const Amg::Vector3D& v2) {
	double dp = v1.dot(v2);
	dp /= v1.mag() * v2.mag();
	if (dp > 1)
		dp = 1;
	if (dp < -1)
		dp = -1;
	return acos(dp);
}


/** calculates the squared distance between two point in 3D space */
inline float distance2(const Amg::Vector3D& p1, const Amg::Vector3D& p2) {
	float dx = p2.x()-p1.x(), dy = p2.y()-p1.y(), dz = p2.z()-p1.z();
	return dx*dx + dy*dy + dz*dz;
}

/** calculates the distance between two point in 3D space */
inline float distance(const Amg::Vector3D& p1, const Amg::Vector3D& p2) {
	return std::sqrt( distance2(p1, p2) );
}




/** sets the phi angle of a vector without changing theta nor the magnitude */
inline void setPhi(Amg::Vector3D& v, double phi) {
	double xy = v.perp();
        v[0] = xy * cos(phi);
        v[1] = xy * sin(phi);
}

/** sets the theta and phi angle of a vector without changing the magnitude */
inline void setThetaPhi(Amg::Vector3D& v, double theta, double phi) {
	double mag = v.mag();
    	v[0] = mag * sin(theta) * cos(phi);
    	v[1] = mag * sin(theta) * sin(phi);
    	v[2] = mag * cos(theta);
}

/** sets radius, the theta and phi angle of a vector. Angles are measured in RADIANS */
inline void setRThetaPhi(Amg::Vector3D& v, double r, double theta, double phi) {
    v[0] = r * sin(theta) * cos(phi);
    v[1] = r * sin(theta) * sin(phi);
    v[2] = r * cos(theta);
}

/** sets the theta of a vector without changing phi nor the magnitude */
inline void setTheta(Amg::Vector3D& v, double theta) {
	setThetaPhi(v, theta, v.phi());
}

/** scales the vector in the xy plane without changing the z coordinate nor the angles */
inline void setPerp(Amg::Vector3D& v, double perp) {
	double p = v.perp();
	if (p != 0.0) {
		double scale = perp / p;
		v[0] *= scale;
		v[1] *= scale;
	}
}

/** scales the vector length without changing the angles */
inline void setMag(Amg::Vector3D& v, double mag) {
	double p = v.mag();
	if (p != 0.0) {
		double scale = mag / p;
		v[0] *= scale;
		v[1] *= scale;
		v[2] *= scale;
	}
}
inline double deltaPhi(const Amg::Vector3D& v1, const Amg::Vector3D& v2) {
	double dphi = v2.phi() - v1.phi();
	if (dphi > M_PI) {
		dphi -= M_PI*2;
	} else if (dphi <= -M_PI) {
		dphi += M_PI*2;
	}
	return dphi;
}
inline double deltaR(const Amg::Vector3D& v1, const Amg::Vector3D& v2){
	double a = v1.eta() - v2.eta();
	double b = deltaPhi(v1,v2);
	return sqrt(a*a + b*b);
}






/**
 * Sets components in cartesian coordinate system.
 */
inline void setVector3DCartesian(Amg::Vector3D& v1, double x1, double y1, double z1) { v1[0] = x1; v1[1] = y1; v1[2] = z1; }
/**
 * Gets magnitude squared of the vector.
 */
inline double mag2Vector3D(const Amg::Vector3D& v1) { return v1.x()*v1.x() + v1.y()*v1.y() + v1.z()*v1.z(); }
/**
 * Gets magnitude of the vector.
 */
inline double magVector3D(const Amg::Vector3D& v1) { return std::sqrt(mag2Vector3D(v1)); }
/**
 * Gets r-component in spherical coordinate system
 */
inline double rVector3D(const Amg::Vector3D& v1) { return magVector3D(v1); }

/**
 * Transform a point from a Trasformation3D
 *
 * from CLHEP::Point3D::transform:
 * http://proj-clhep.web.cern.ch/proj-clhep/doc/CLHEP_2_0_4_7/doxygen/html/Point3D_8cc-source.html#l00032
 */
inline Amg::Vector3D transform( Amg::Vector3D& v, Amg::Transform3D& tr ) {
	Amg::Vector3D vect;
	double vx = v.x(), vy = v.y(), vz = v.z();
	setVector3DCartesian( vect,
			tr(0,0)*vx + tr(0,1)*vy + tr(0,2)*vz + tr(0,3),
			tr(1,0)*vx + tr(1,1)*vy + tr(1,2)*vz + tr(1,3),
			tr(2,0)*vx + tr(2,1)*vy + tr(2,2)*vz + tr(2,3));
	return vect;
}




/*
 * the analogous to CLHEP HepGeom::Transform3D trans (localRot, theSurface.transform().translation());
 */
inline Amg::Transform3D getTransformFromRotTransl(Amg::RotationMatrix3D rot, Amg::Vector3D transl_vec )
{
	Amg::Transform3D trans = Amg::Transform3D::Identity();
    trans = trans * rot;
    trans.translation() = transl_vec;
	return trans;
}

/*
 * Replacing the CLHEP::HepRotation::getAngleAxis() functionality
 *
 * Note:
 * CLHEP has a 'HepRotation::getAngleAxis()' function, e.g.:
 * ---
 * CLHEP::HepRotation rotation   = transform.getRotation();
 * CLHEP::Hep3Vector  rotationAxis;
 * double      rotationAngle;
 * rotation.getAngleAxis(rotationAngle,rotationAxis);
 * ---
 */
inline void getAngleAxisFromRotation(Amg::RotationMatrix3D& rotation, double& rotationAngle, Amg::Vector3D& rotationAxis)
{
	rotationAngle = 0.;

	double xx = rotation(0,0);
	double yy = rotation(1,1);
	double zz = rotation(2,2);

	double cosa  = 0.5 * (xx + yy + zz - 1);
	double cosa1 = 1 - cosa;

	if (cosa1 <= 0) {
		rotationAngle = 0;
		rotationAxis  = Amg::Vector3D(0,0,1);
	}
	else{
		double x=0, y=0, z=0;
		if (xx > cosa) x = sqrt((xx-cosa)/cosa1);
		if (yy > cosa) y = sqrt((yy-cosa)/cosa1);
		if (zz > cosa) z = sqrt((zz-cosa)/cosa1);
		if (rotation(2,1) < rotation(1,2)) x = -x;
		if (rotation(0,2) < rotation(2,0)) y = -y;
		if (rotation(1,0) < rotation(0,1)) z = -z;
		rotationAngle = (cosa < -1.) ? acos(-1.) : acos(cosa);
		rotationAxis  = Amg::Vector3D(x,y,z);
	}

	return;
}

/**
 * Get the Translation vector out of a Transformation
 */
inline Amg::Vector3D getTranslationVectorFromTransform(const Amg::Transform3D& tr) {
	return Amg::Vector3D(tr(0,3),tr(1,3),tr(2,3));
} // TODO: check! it's perhaps useless, you acn use the transform.translation() method



/**
 * get a AngleAxis from an angle and an axis.
 *
 * to replace the CLHEP constructor:
 * CLHEP::Rotate3D::Rotate3D(double a, cconst Vector3D< double > & v)
 */
inline Amg::Rotation3D getRotation3DfromAngleAxis(double angle, Amg::Vector3D& axis)
{
	AngleAxis3D t;
	t = Eigen::AngleAxis<double>(angle,axis);

	Amg::Rotation3D rot;
	rot = t;

	return rot;
}


/**
 * get a rotation transformation around X-axis
 */
inline Amg::Transform3D getRotateX3D(double angle) {
	Amg::Transform3D transf;
	Amg::AngleAxis3D angleaxis(angle, Amg::Vector3D(1.,0.,0.));
	transf = angleaxis;
	return transf;
}
/**
 * get a rotation transformation around Y-axis
 */
inline Amg::Transform3D getRotateY3D(double angle) {
	Amg::Transform3D transf;
	Amg::AngleAxis3D angleaxis(angle, Amg::Vector3D(0.,1.,0.));
	transf = angleaxis;
	return transf;
}
/**
 * get a rotation transformation around Z-axis
 */
inline Amg::Transform3D getRotateZ3D(double angle) {
	Amg::Transform3D transf;
	Amg::AngleAxis3D angleaxis(angle, Amg::Vector3D(0.,0.,1.));
	transf = angleaxis;
	return transf;
}



} // end of Amg namespace

#endif
