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
Acts::ExtrapolationCell<T>::stepTransport(std::unique_ptr<const T>                 stepParameters,
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
  if (stepConfig.checkMode(ExtrapolationMode::Destination)){
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
Acts::ExtrapolationCell<T>::stepMaterial(std::unique_ptr<const T>  stepParameters,
                                         const Vector3D&           stepPosition,
                                         const Surface&            stepSurface,
                                         double                    stepFactor,
                                         const MaterialProperties* mprop)
{
  // if there's new stepParameters then change the lead
  if (stepParameters){
    lastLeadParameters = leadParameters;
    leadParameters     = stepParameters.get();
  }
  // add material to the global counters
  if (mprop) {
    // the overal material
    materialX0 += stepFactor * mprop->thicknessInX0();
    materialL0 += stepFactor * mprop->thicknessInL0();
  }  
  // prepare the extrapolation mode
  std::vector<ExtrapolationMode::eMode> emode = {ExtrapolationMode::CollectMaterial};
  // check the last step if this is a potential update on
  if (stepParameters &&
      extrapolationSteps.size() && 
      extrapolationSteps.back().configuration.checkMode(ExtrapolationMode::Destination)){
      // we move the endParameters into the last step
      // this should set endParameters to nullptr
      extrapolationSteps.back().parameters = std::move(endParameters);
      // set the new endParameters to the stepParameters
      endParameters = std::move(stepParameters);
      // update the current mode to also be Destination
      emode.push_back(ExtrapolationMode::Destination);
  }
  ExtrapolationConfig stepConfig(emode);
  extrapolationSteps.push_back(ExtrapolationStep<T>(std::move(stepParameters),
                                                    &stepSurface,
                                                    stepConfig,
                                                    mprop,
                                                    nullptr));
  // complete the step information
  extrapolationSteps.back().position         = stepPosition;                                                  
  extrapolationSteps.back().materialScaling  = stepFactor;                                                    
                                                    

}

