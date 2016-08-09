// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericLayerBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Utilities/Helpers.hpp"
#include "ACTS/Detector/DetectorElementBase.hpp"
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Utilities/ApproachDescriptor.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinnedArray1D.hpp"
#include "ACTS/Utilities/BinnedArray2D.hpp"
#include "ACTS/Examples/GenericDetectorElement.hpp"
#include "ACTS/Examples/GenericLayerBuilder.hpp"
#include <iostream>

Acts::GenericLayerBuilder::GenericLayerBuilder(
    const Acts::GenericLayerBuilder::Config& glbConfig)
  : m_nLayers(), m_cLayers(), m_pLayers()
{
  ACTS_DEBUG("initialize()");
  /// @TODO a configuraiton check should be done here
  setConfiguration(glbConfig);
  // Tool needs to be initialized
  constructLayers();
}

void
Acts::GenericLayerBuilder::setConfiguration(
    const Acts::GenericLayerBuilder::Config& glbConfig)
{
  // @TODO check consistency
  // copy the configuration
  m_cfg = glbConfig;
}

Acts::GenericLayerBuilder::~GenericLayerBuilder()
{
}

void
Acts::GenericLayerBuilder::constructLayers()
{
  
  size_t imodule = 0;
  // ----------------------- central layers -------------------------
  // the central layers
  size_t numcLayers = m_cfg.centralLayerRadii.size();
  if (numcLayers) {
    ACTS_DEBUG("Configured to build " << numcLayers << " active central layers.");
    m_cLayers.reserve(numcLayers);
    // loop through
    for (size_t icl = 0; icl < numcLayers; ++icl) {
      // layer R/Z
      double layerR = m_cfg.centralLayerRadii.at(icl);
      // some screen output
      ACTS_VERBOSE("Build layer " << icl << " with target radius = " << layerR);
      
      // prepare the Surface vector
      std::vector<const Surface*> sVector;
      // assign the current envelope
      double layerEnvelopeCoverZ = m_cfg.centralLayerEnvelopes.size()
          ? m_cfg.centralLayerEnvelopes.at(icl).second : 0.;
      // module size & tilt
      double modulePhiTilt   = m_cfg.centralModuleTiltPhi.at(icl);
      double moduleHalfX     = m_cfg.centralModuleHalfX.at(icl);
      double moduleHalfY     = m_cfg.centralModuleHalfY.at(icl);
      double moduleThickness = m_cfg.centralModuleThickness.at(icl);
      // create the shared module
      std::shared_ptr<const PlanarBounds> moduleBounds(
          new RectangleBounds(moduleHalfX, moduleHalfY));
      // Identifier @TODO unique Identifier - use a GenericDetector identifier
      size_t nCetralModules = 
          m_cfg.centralModuleBinningSchema.at(icl).first
        * m_cfg.centralModuleBinningSchema.at(icl).second;
      
      ACTS_VERBOSE("- number of modules " << nCetralModules 
       << " ( from " << m_cfg.centralModuleBinningSchema.at(icl).first << " x "
       << m_cfg.centralModuleBinningSchema.at(icl).second << " )");
      
      sVector.reserve(nCetralModules);
      // create the Module material from input
      std::shared_ptr<const SurfaceMaterial> moduleMaterialPtr = nullptr;
      if (m_cfg.centralModuleMaterial.size()) {
        // get the sensor material from configuration
        Material moduleMaterial = m_cfg.centralModuleMaterial.at(icl);
        MaterialProperties moduleMaterialProperties(moduleMaterial,
                                                    moduleThickness);
        // create a new surface material                                            
        moduleMaterialPtr = std::shared_ptr<const SurfaceMaterial>(
            new HomogeneousSurfaceMaterial(moduleMaterialProperties));
      }

      // confirm
      if (m_cfg.centralModulePositions.at(icl).size() != nCetralModules){
         ACTS_WARNING("Mismatching module numbers, configuration error!");
         ACTS_WARNING("- Binning schema suggests : " << nCetralModules);
         ACTS_WARNING("- Positions provided are  : " << m_cfg.centralModulePositions.at(icl).size());
         
       }
      // loop over the position, create the modules 
      for (auto& moduleCenter : m_cfg.centralModulePositions.at(icl)){
        // create the association transform
        double modulePhi = moduleCenter.phi();
        // the local z axis is the normal vector
        Vector3D moduleLocalZ(cos(modulePhi + modulePhiTilt),
                              sin(modulePhi + modulePhiTilt),
                              0.);
        // the local y axis is the global z axis                      
        Vector3D moduleLocalY(0., 0., 1);
        // the local x axis the normal to local y,z
        Vector3D moduleLocalX(-sin(modulePhi + modulePhiTilt),
                              cos(modulePhi + modulePhiTilt),
                              0.);
        // create the RotationMatrix
        RotationMatrix3D moduleRotation;
        moduleRotation.col(0) = moduleLocalX;
        moduleRotation.col(1) = moduleLocalY;
        moduleRotation.col(2) = moduleLocalZ;
        // get the moduleTransform
        std::shared_ptr<Transform3D> moduleTransform(new Transform3D(
            getTransformFromRotTransl(moduleRotation, moduleCenter)));
        // stereo angle if necessary
        if (m_cfg.centralModuleFrontsideStereo.size()
            && m_cfg.centralModuleFrontsideStereo.at(icl) != 0.) {
          // twist by the stereo angle
          double stereo = m_cfg.centralModuleFrontsideStereo.at(icl);
          (*moduleTransform.get()) *= AngleAxis3D(-stereo, Vector3D::UnitZ());
        }
        // count the modules
        ++imodule;
        Identifier moduleIdentifier
            = Identifier(Identifier::value_type(imodule));
        // create the module
        DetectorElementBase* module
            = new GenericDetectorElement(moduleIdentifier,
                                         moduleTransform,
                                         moduleBounds,
                                         moduleThickness,
                                         moduleMaterialPtr);
        // register the surface
        sVector.push_back(&module->surface());                              
        // store the module
        // @TODO detector store facility
        m_centralModule.push_back(module);
        // IF double modules exist
        // and the backside one (if configured to do so)
        if (m_cfg.centralModuleBacksideGap.size()) {
          // ncrease the counter @TODO switch to identifier service
          ++imodule;
          // create the module identifier
          moduleIdentifier = Identifier(Identifier::value_type(imodule));
          moduleCenter     = moduleCenter
              + m_cfg.centralModuleBacksideGap.at(icl) * moduleLocalZ;
          moduleTransform = std::shared_ptr<Transform3D>(new Transform3D(
              getTransformFromRotTransl(moduleRotation, moduleCenter)));
          // apply the stereo
          if (m_cfg.centralModuleBacksideStereo.size()) {
            // twist by the stereo angle
            double stereoBackSide
                = m_cfg.centralModuleBacksideStereo.at(icl);
            (*moduleTransform.get())
                *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
          }
          // everything is set for the next module
          DetectorElementBase* bsmodule
              = new GenericDetectorElement(moduleIdentifier,
                                           moduleTransform,
                                           moduleBounds,
                                           moduleThickness,
                                           moduleMaterialPtr);
          // register the backside as bin member
          std::vector<const DetectorElementBase*> bsbinmember = {module};
          std::vector<const DetectorElementBase*> binmember   = {bsmodule};
          bsmodule->registerBinmembers(bsbinmember);
          module->registerBinmembers(binmember);
          // memory management - we need a detector store to hold them
          // somewhere @TODO detector store facility
          m_centralModule.push_back(bsmodule);
        }
      }
      
      size_t phiBins  = m_cfg.centralModuleBinningSchema.at(icl).first;
      phiBins *= m_cfg.centralLayerBinMultipliers.first;
      size_t zBins    = m_cfg.centralModuleBinningSchema.at(icl).second;
      zBins *= m_cfg.centralLayerBinMultipliers.second;
      // create the surface array - it will also fill the accesible binmember
      // chache if avalable
      LayerPtr cLayer = m_cfg.layerCreator->cylinderLayer(
          sVector,
          m_cfg.approachSurfaceEnvelope,
          layerEnvelopeCoverZ,
          phiBins,
          zBins);
      // the layer is built le't see if it needs material
      if (m_cfg.centralLayerMaterialProperties.size()) {
        // get the material from configuration
        MaterialProperties layerMaterialProperties =
          m_cfg.centralLayerMaterialProperties.at(icl);
        std::shared_ptr<const SurfaceMaterial> layerMaterialPtr(
            new HomogeneousSurfaceMaterial(layerMaterialProperties));
        // central material
        if (m_cfg.centralLayerMaterialConcentration.at(icl) == 0.) {
          // the layer surface is the material surface
          cLayer->surfaceRepresentation().setAssociatedMaterial(layerMaterialPtr);
        } else {
          // approach surface material
          // get the approach descriptor - at this stage we know that the
          // approachDescriptor exists
          auto approachSurfaces
              = cLayer->approachDescriptor()->containedSurfaces();
          if (m_cfg.centralLayerMaterialConcentration.at(icl) > 0) {
            approachSurfaces.at(1)->setAssociatedMaterial(layerMaterialPtr);
              ACTS_VERBOSE("- and material at outer approach surface");
          } else {
            approachSurfaces.at(0)->setAssociatedMaterial(layerMaterialPtr);
            ACTS_VERBOSE("- and material at inner approach surface");
          }
        }
      }
      // push it into the layer vector
      m_cLayers.push_back(cLayer);
    }
  }


  // -------------------------------- endcap type layers
  // pos/neg layers
  size_t numpnLayers = m_cfg.posnegLayerPositionsZ.size();
  if (numpnLayers) {
    ACTS_DEBUG("Configured to build 2 * " << numpnLayers
               << " passive positive/negative side layers.");
    m_pLayers.reserve(numpnLayers);
    m_nLayers.reserve(numpnLayers);
    
    for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl) {
      // some screen output
      ACTS_VERBOSE("- build layers " << (2*ipnl) << " and "
        <<  (2*ipnl)+1 << "at +/- z = " << m_cfg.posnegLayerPositionsZ.at(ipnl));    
      /// some preparation work
      // define the layer envelope
      double layerEnvelopeR = m_cfg.posnegLayerEnvelopeR.at(ipnl);
      // prepare for the r binning
      std::vector<const Surface*> nsVector;
      std::vector<const Surface*> psVector;
      // now fill the vectors 
      for (auto& discModulePositions : m_cfg.posnegModulePositions.at(ipnl)){
        size_t ipnR = 0;
        for (auto& ringModulePosition : discModulePositions){
          // module specifications
          double moduleThickness = m_cfg.posnegModuleThickness.at(ipnl).at(ipnR);
          double moduleMinHalfX  = m_cfg.posnegModuleMinHalfX.at(ipnl).at(ipnR);
          double moduleMaxHalfX  = m_cfg.posnegModuleMaxHalfX.size() ?
          m_cfg.posnegModuleMaxHalfX.at(ipnl).at(ipnR) : 0.;
          double moduleHalfY     = m_cfg.posnegModuleHalfY.at(ipnl).at(ipnR);
          // module material
          // create the Module material from input
          std::shared_ptr<const SurfaceMaterial> moduleMaterialPtr = nullptr;
          if (m_cfg.posnegModuleMaterial.size()) {
            MaterialProperties moduleMaterialProperties(m_cfg.posnegModuleMaterial.at(ipnl).at(ipnR),
                                                        moduleThickness);
            // and create the shared pointer
            moduleMaterialPtr = std::shared_ptr<const SurfaceMaterial>(
                new HomogeneousSurfaceMaterial(moduleMaterialProperties));
          } 
          // create the bounds
          PlanarBounds* pBounds = nullptr;
          if (moduleMaxHalfX != 0. && moduleMinHalfX != moduleMaxHalfX)
            pBounds = new TrapezoidBounds(
                moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
          else
            pBounds = new RectangleBounds(moduleMaxHalfX, moduleHalfY);
          // now create the shared bounds from it
          std::shared_ptr<const PlanarBounds> moduleBounds(pBounds);
          // the module transform from the position
          double modulePhi = ringModulePosition.phi();
          // the center position of the modules
          Vector3D pModuleCenter(ringModulePosition);
          // take the mirrored position wrt x/y
          Vector3D nModuleCenter(pModuleCenter.x(),
                                 pModuleCenter.y(),
                                 -pModuleCenter.z());
          // the rotation matrix of the module
          Vector3D moduleLocalY(cos(modulePhi), sin(modulePhi), 0.);
          // take different axis to have the same readout direction
          Vector3D pModuleLocalZ(0., 0., 1.);  
          // take different axis to have the same readout direction
          Vector3D nModuleLocalZ(0., 0., -1.);  
          Vector3D nModuleLocalX = moduleLocalY.cross(nModuleLocalZ);
          Vector3D pModuleLocalX = moduleLocalY.cross(pModuleLocalZ);
          // local rotation matrices
          // create the RotationMatrix - negative side
          RotationMatrix3D nModuleRotation;
          nModuleRotation.col(0) = nModuleLocalX;
          nModuleRotation.col(1) = moduleLocalY;
          nModuleRotation.col(2) = nModuleLocalZ;
          // create the RotationMatrix - positive side
          RotationMatrix3D pModuleRotation;
          pModuleRotation.col(0) = pModuleLocalX;
          pModuleRotation.col(1) = moduleLocalY;
          pModuleRotation.col(2) = pModuleLocalZ;
          // the transforms for the two modules
          std::shared_ptr<Transform3D> nModuleTransform(new Transform3D(
              getTransformFromRotTransl(nModuleRotation, nModuleCenter)));
          std::shared_ptr<Transform3D> pModuleTransform(new Transform3D(
              getTransformFromRotTransl(pModuleRotation, pModuleCenter)));
          // create the modules identifier @TODO Idenfier service
          Identifier nModuleIdentifier
              = Identifier(Identifier::value_type(2 * imodule));
          Identifier pModuleIdentifier
              = Identifier(Identifier::value_type(2 * imodule + 1));
          // create the module
          GenericDetectorElement* nmodule
              = new GenericDetectorElement(nModuleIdentifier,
                                           nModuleTransform,
                                           moduleBounds,
                                           moduleThickness,
                                           moduleMaterialPtr);
          GenericDetectorElement* pmodule
              = new GenericDetectorElement(pModuleIdentifier,
                                           pModuleTransform,
                                           moduleBounds,
                                           moduleThickness,
                                           moduleMaterialPtr);
          // memory management - we need a detector store to hold them somewhere
          // @TODO add detector store facility
          m_posnegModule.push_back(nmodule);
          m_posnegModule.push_back(pmodule);
          // now deal with the potential backside
          if (m_cfg.posnegModuleBacksideGap.size()) {
            // ncrease the counter @TODO switch to identifier service
            nModuleIdentifier = Identifier(Identifier::value_type(++imodule));
            pModuleIdentifier = Identifier(Identifier::value_type(++imodule));
            // the new centers
            nModuleCenter = nModuleCenter
                + m_cfg.posnegModuleBacksideGap.at(ipnl).at(ipnR)
                    * nModuleLocalZ;
            pModuleCenter = pModuleCenter
                + m_cfg.posnegModuleBacksideGap.at(ipnl).at(ipnR)
                    * pModuleLocalZ;
            // the new transforms
            nModuleTransform = std::shared_ptr<Transform3D>(new Transform3D(
                getTransformFromRotTransl(nModuleRotation, nModuleCenter)));
            pModuleTransform = std::shared_ptr<Transform3D>(new Transform3D(
                getTransformFromRotTransl(pModuleRotation, pModuleCenter)));
            // apply the stereo
            if (m_cfg.posnegModuleBacksideStereo.size()) {
              // twist by the stereo angle
              double stereoBackSide
                  = m_cfg.posnegModuleBacksideStereo.at(ipnl).at(ipnR);
              (*nModuleTransform.get())
                  *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
              (*pModuleTransform.get())
                  *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
            }
            // everything is set for the next module
            GenericDetectorElement* bsnmodule
                = new GenericDetectorElement(nModuleIdentifier,
                                             nModuleTransform,
                                             moduleBounds,
                                             moduleThickness,
                                             moduleMaterialPtr);
            GenericDetectorElement* bspmodule
                = new GenericDetectorElement(pModuleIdentifier,
                                             pModuleTransform,
                                             moduleBounds,
                                             moduleThickness,
                                             moduleMaterialPtr);
            // register the backside of the binmembers
            std::vector<const DetectorElementBase*> bspbinmember = {pmodule};
            std::vector<const DetectorElementBase*> pbinmember   = {bspmodule};
            std::vector<const DetectorElementBase*> bsnbinmember = {nmodule};
            std::vector<const DetectorElementBase*> nbinmember   = {bsnmodule};
            bsnmodule->registerBinmembers(bsnbinmember);
            nmodule->registerBinmembers(nbinmember);
            bspmodule->registerBinmembers(bspbinmember);
            pmodule->registerBinmembers(pbinmember);
            // memory management - we need a detector store to hold them
            // somewhere @TODO add detector store facility
            m_posnegModule.push_back(bsnmodule);
            m_posnegModule.push_back(bspmodule);
          }
          // create the surface
          nsVector.push_back(&nmodule->surface());
          psVector.push_back(&pmodule->surface());
        }
        // counter of rings
        ++ipnR;
      }
      // @TODO this needs an update
      size_t layerBinsR   = m_cfg.posnegModulePhiBins.at(ipnl).size();
      size_t layerBinsPhi = m_cfg.posnegModulePhiBins.at(ipnl).at(0);
      // create the layers with the surface arrays
      LayerPtr nLayer
          = m_cfg.layerCreator->discLayer(nsVector,
                                          layerEnvelopeR,
                                          layerEnvelopeR,
                                          m_cfg.approachSurfaceEnvelope,
                                          layerBinsR,
                                          layerBinsPhi);
      LayerPtr pLayer
          = m_cfg.layerCreator->discLayer(psVector,
                                          layerEnvelopeR,
                                          layerEnvelopeR,
                                          m_cfg.approachSurfaceEnvelope,
                                          layerBinsR,
                                          layerBinsPhi);
                                          
      // the layer is built le't see if it needs material
      if (m_cfg.posnegLayerMaterialProperties.size()) {
          std::shared_ptr<const SurfaceMaterial> layerMaterialPtr(
            new HomogeneousSurfaceMaterial(m_cfg.posnegLayerMaterialProperties[ipnl]));
        // central material
        if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) == 0.) {
          // assign the surface material - the layer surface is the material
          // surface
          nLayer->surfaceRepresentation().setAssociatedMaterial(layerMaterialPtr);
          pLayer->surfaceRepresentation().setAssociatedMaterial(layerMaterialPtr);
        } else {
          // approach surface material
          // get the approach descriptor - at this stage we know that the
          // approachDescriptor exists
          auto nApproachSurfaces
              = nLayer->approachDescriptor()->containedSurfaces();
          auto pApproachSurfaces
              = pLayer->approachDescriptor()->containedSurfaces();
          if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) > 0.) {
            nApproachSurfaces.at(0)->setAssociatedMaterial(layerMaterialPtr);
            pApproachSurfaces.at(1)->setAssociatedMaterial(layerMaterialPtr);
          } else {
            nApproachSurfaces.at(1)->setAssociatedMaterial(layerMaterialPtr);
            pApproachSurfaces.at(0)->setAssociatedMaterial(layerMaterialPtr);
          }
        }
      }
      // push it into the layer vector
      m_nLayers.push_back(nLayer);
      m_pLayers.push_back(pLayer);                                          
    }
  }
}
