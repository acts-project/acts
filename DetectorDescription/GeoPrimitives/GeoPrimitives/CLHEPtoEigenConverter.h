/*
 * CLHEPtoEigenConverter.h
 *
 *  Created on: Jun 26, 2013
 *      Author: rlangenb
 *
 *      update: rbianchi - Feb 25 2014
 */

#ifndef CLHEPTOEIGENCONVERTER_H_
#define CLHEPTOEIGENCONVERTER_H_

// Make it clear that this header is not for standalone usage:
#ifdef XAOD_STANDALONE
#error "This header is not meant to be used in standalone mode"
#endif // XAOD_STANDALONE

#include "GeoPrimitives/GeoPrimitives.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace Amg {

    inline Amg::Transform3D CLHEPTransformToEigen(
            const HepGeom::Transform3D& CLHEPtransf) {
        Amg::Transform3D t = Amg::Transform3D();
        //loop unrolled for performance
        t(0, 0) = CLHEPtransf(0, 0);
        t(0, 1) = CLHEPtransf(0, 1);
        t(0, 2) = CLHEPtransf(0, 2);
        t(1, 0) = CLHEPtransf(1, 0);
        t(1, 1) = CLHEPtransf(1, 1);
        t(1, 2) = CLHEPtransf(1, 2);
        t(2, 0) = CLHEPtransf(2, 0);
        t(2, 1) = CLHEPtransf(2, 1);
        t(2, 2) = CLHEPtransf(2, 2);
        t(0, 3) = CLHEPtransf(0, 3);
        t(1, 3) = CLHEPtransf(1, 3);
        t(2, 3) = CLHEPtransf(2, 3);
        return t;
    }
    
    inline RotationMatrix3D CLHEPRotationToEigen(
            const CLHEP::HepRotation& CLHEProtation) {
        Amg::RotationMatrix3D t;

        t(0, 0) = CLHEProtation(0, 0);
        t(0, 1) = CLHEProtation(0, 1);
        t(0, 2) = CLHEProtation(0, 2);
        t(1, 0) = CLHEProtation(1, 0);
        t(1, 1) = CLHEProtation(1, 1);
        t(1, 2) = CLHEProtation(1, 2);
        t(2, 0) = CLHEProtation(2, 0);
        t(2, 1) = CLHEProtation(2, 1);
        t(2, 2) = CLHEProtation(2, 2);
        return t;
    }
    inline Translation3D CLHEPTranslationToEigen(
            const CLHEP::Hep3Vector& CLHEPtranslation) {
        return Translation3D(
                Vector3D(CLHEPtranslation[0], CLHEPtranslation[1],
                        CLHEPtranslation[2]));
    }


    // from: http://proj-clhep.web.cern.ch/proj-clhep/doc/CLHEP_2_0_4_7/doxygen/html/classHepGeom_1_1Translate3D.html#f2df65781931c7df9cc2858de2c89151
	//Amg::Transform3D(1, 0, 0, CLHEPtranslate3D[0],
	//                 0, 1, 0, CLHEPtranslate3D[1],
	//  		       0, 0, 1, CLHEPtranslate3D[2]);
    inline Amg::Transform3D CLHEPTranslate3DToEigen(
            const HepGeom::Translate3D& CLHEPtranslate3D)
    {
    	Amg::Transform3D t = Amg::Transform3D();
    	t.setIdentity();
    	t(0, 3) = CLHEPtranslate3D(0, 3);
    	t(1, 3) = CLHEPtranslate3D(1, 3);
    	t(2, 3) = CLHEPtranslate3D(2, 3);
        return t;
    }
    
    inline HepGeom::Transform3D EigenTransformToCLHEP(
            const Amg::Transform3D& eigenTransf) {
        CLHEP::HepRotation rotation(
                CLHEP::Hep3Vector(eigenTransf(0, 0), eigenTransf(1, 0), eigenTransf(2, 0)),
                CLHEP::Hep3Vector(eigenTransf(0, 1), eigenTransf(1, 1), eigenTransf(2, 1)),
                CLHEP::Hep3Vector(eigenTransf(0, 2), eigenTransf(1, 2), eigenTransf(2, 2)));
        CLHEP::Hep3Vector translation(eigenTransf(0, 3), eigenTransf(1, 3), eigenTransf(2, 3));
        HepGeom::Transform3D t(rotation, translation);
        return t;
    }

    inline Amg::Vector3D Hep3VectorToEigen(const CLHEP::Hep3Vector& CLHEPvector) {
        return Vector3D(CLHEPvector[0], CLHEPvector[1], CLHEPvector[2]);
    }

    inline CLHEP::Hep3Vector EigenToHep3Vector(const Amg::Vector3D& eigenvector) {
        return CLHEP::Hep3Vector(eigenvector[0], eigenvector[1], eigenvector[2]);
    }
}

#endif
