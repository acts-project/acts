// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationCell.ipp, ACTS project
///////////////////////////////////////////////////////////////////

template <class T>
void
Acts::ExtrapolationCell<T>::stepTransport(
    std::unique_ptr<const T>                 stepParameters,
    std::vector<ExtrapolationMode::eMode>    stepModes,
    double                                   stepLength,
    std::unique_ptr<const TransportJacobian> stepJacobian)
{
  // remember the last lead parameters
  lastLeadParameters = leadParameters;
  // add the path length to the global counter
  pathLength += stepLength;
  // these are new parameters created by transport/propagation
  // set them as new lead parameters
  leadParameters = stepParameters.get();
  // current step surface
  const Surface* stepSurface = &(stepParameters->referenceSurface());
  // create a configuration for this step holding all modes
  ExtrapolationConfig stepConfig = ExtrapolationConfig(stepModes);
  // check if we have the destination
  if (stepConfig.checkMode(ExtrapolationMode::Destination)) {
    // this should set the stepParameters to nullptr
    endParameters = std::move(stepParameters);
  }
  // this is a new step with , so fill it
  extrapolationSteps.push_back(ExtrapolationStep<T>(std::move(stepParameters),
                                                    stepSurface,
                                                    stepConfig,
                                                    nullptr,
                                                    std::move(stepJacobian),
                                                    stepLength));
}

template <class T>
void
Acts::ExtrapolationCell<T>::stepMaterial(
    std::unique_ptr<const T>  stepParameters,
    const Vector3D&           stepPosition,
    const Surface&            stepSurface,
    double                    stepFactor,
    const MaterialProperties* mprop)
{
  
  // if this is on a new surface, 
  // so create a new extrapolation step
  if (extrapolationSteps.size() && 
      (&stepSurface) != extrapolationSteps.back().surface){
      // if the last step was already the destination, 
      // then do not create a final step
      // let's check the last one 
      auto& lstep = extrapolationSteps.back();
      if (lstep.configuration.checkMode(ExtrapolationMode::Destination)){
        // in case the last propagation was towards a Desitination
        // then they are written into the endParameters
        // after additional material update, we move the endParameters
        // into the last step, which sets endParameters to nullptr
        extrapolationSteps.back().parameters = std::move(endParameters);
        // then set the new endParameters to the be the stepParameters
        // this sets the stepParameters to nullptr
        endParameters = std::move(stepParameters);
      } else {
        // create a new destination  
        extrapolationSteps.push_back(ExtrapolationStep<T>());
      }
  } else {
    // a new step is needed for the first one in any case
    extrapolationSteps.push_back(ExtrapolationStep<T>());
  }

  // we work with the last one it's either 
  // - a nelwy created one
  // - the last one 
  auto& cstep = extrapolationSteps.back();

  // if there's new stepParameters then change the lead
  if (stepParameters) {
    // store the old parameters if present
    if (cstep.parameters) 
        cstep.preparameters = std::move(cstep.parameters);
    // bookkeeping    
    lastLeadParameters = leadParameters;
    leadParameters     = stepParameters.get();
    // setting
    cstep.parameters   = std::move(stepParameters);
  }
  // add material to the global counters
  if (mprop) {
    // the overal material
    materialX0 += stepFactor * mprop->thicknessInX0();
    materialL0 += stepFactor * mprop->thicknessInL0();
  }
  // simply add the material configuration
  cstep.configuration.addMode(ExtrapolationMode::CollectMaterial);
  // finalise the step information
  cstep.surface         = &stepSurface;
  cstep.material        = mprop;
  cstep.position        = stepPosition;
  cstep.materialScaling = stepFactor;

}
