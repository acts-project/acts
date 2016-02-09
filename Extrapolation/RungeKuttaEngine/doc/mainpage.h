/**
@mainpage The TrkExRungeKutta package

Propagation of TrackParameters and their associated covariance matrices in a Runge-Kutta integrated method

@section TrkExRkIntroduction Introduction

The RungeKuttaEngine implements the IPropagator interface with a Runge-Kutta 
integrated method to be used in a non-homogenous magnetic field.
In case of absence of a magnetic field the RungeKuttaEngine would propagate the track-parameters
according to a straight line, in case of an homogenous magnetic field, a helical propagation is performed.

Common transformations with the STEP_Propagator can be found in the ExtrapolationUtils package.

The tool implements two abstract interfaces:
- under Ats::IPropagator it provides propagation of EDM Ats::TrackParameters
- under Ats::IPatternParametersPropagator it provides propagation of parameters for internal use, the Ats::PatternTrackParameters

@section TrkExRkComments Comments

Please let me know of any errors, or if anything is unclear.
@author Igor.Gavrilenko@cern.ch

*/

