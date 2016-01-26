/*
 * eigen_migration
 * AmgStringHelpers.h
 *
 *  Created on: Feb 25, 2014
 *      Author: rbianchi <Riccardo.Maria.Bianchi@cern.ch>
 *
 */

#ifndef AMGSTRINGHELPERS_H_
#define AMGSTRINGHELPERS_H_

#include "GeoPrimitives/GeoPrimitives.h"

namespace Amg{


/**
 * write an Amg Eigen object to std::string
 */
template<class T>
std::string AsString(const T& m) {
	std::stringstream tmp;
	tmp << m; // we just use eigen classes capability to write to std::ostream
	return tmp.str();
}



} // end of namespace Amg

#endif /* AMGSTRINGHELPERS_H_ */
