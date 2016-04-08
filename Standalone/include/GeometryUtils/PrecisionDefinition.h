///////////////////////////////////////////////////////////////////
// PrecisionDefinition.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYPRECISION_H
#define ACTS_GEOMETRYUTILS_GEOMETRYPRECISION_H 1

#ifdef TRKDETDESCR_USEFLOATPRECISON 
typedef float TDD_real_t;
#else 
typedef double TDD_real_t;
#endif

#define TDD_max_bound_value 10e10

namespace Acts {
    /** Tolerance for being on Surface */
    static const double s_onSurfaceTolerance = 10e-5;
}

#endif // ACTS_GEOMETRYUTILS_GEOMETRYPRECISION_H
