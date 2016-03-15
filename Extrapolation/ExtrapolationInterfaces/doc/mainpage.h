/**
@mainpage ExtrapolationInterfaces Package

This package contains all interface classes
to be used in the TrkExtrapolation container package.
Concrete implementations can be found in the various propagator packages or in the common
TrkExTools package.

@author Andreas.Salzburger@cern.ch

@section ExtrapolationInterfacesIntro List of interfaces


The following interface classes are defined

\htmlonly
<!-- disable linking: doxygen puts the classes' @brief statement as descriptions to the links and causes mess if that statement has a wikiword inside -->
<ul>
 <li> <b>Acts::IExtrapolator</b>              : interface for the standard Acts::Extrapolator</li>
 <li> <b>Acts::IPropagator</b>                : interface for mathematical propagators</li>
 <li> <b>Acts::IPatternParametersPropagator</b> : interface for propagation of track states given as Acts::PatternTrackParameters</li>
 <li> <b>Acts::IIntersector</b>               : interface for Intersector AlgTools as they are used within iPatRec</li>
 <li> <b>Acts::INavigator</b>                 : interface for the navigator used in the extrapolation process</li>
 <li> <b>Acts::IMaterialEffectsUpdator</b>    : interface for the standard material effects updator</li>
 <li> <b>Acts::IEnergyLossUpdator</b>         : interface for material effects updator (only eloss)</li>
 <li> <b>Acts::IEnergyLossCalculator</b>      : interface for calculating an energy loss correction from a variety of effects</li>
 <li> <b>Acts::IMultipleScatteringUpdator</b> : interface for material effects updator (only multiple scattering)</li>
</ul>
This package is not built as a library, it is a simple include package.   
\endhtmlonly
      
@section ExtrasExtrapolationInterfaces Extra Pages

 - @ref UsedExtrapolationInterfaces
 - @ref reqsExtrapolationInterfaces
*/

/**
@page UsedExtrapolationInterfaces Used Packages
@htmlinclude used_packages.html
*/

/**
@page reqsExtrapolationInterfaces Requirements
@include requirements
*/
