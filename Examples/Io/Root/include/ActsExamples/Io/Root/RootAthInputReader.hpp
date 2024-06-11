// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <mutex>
#include <map>
#include <string>
#include <vector>
#include <memory>

#include "TBranch.h"

class TChain;
namespace ActsExamples{

/// @class RootAthInputReader
///
/// @brief Reader for measurements and spacepoints from an Athena
///        object dumper. 
///        Specifically written for the input ntuple for GNN 
///        See: https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerDetector/InDetGNNTracking/src/DumpObjects.cxx




  
  class RootAthInputReader : public IReader {
  public:
    /// @brief The nested configuration struct
    struct Config {
      // Name of tree
      std::string treename;
      // Name of inputfile
      std::string inputfile;
      // name of the output measurements 
      std::string outputMeasurements = "ath_meas";
      // name of the output pixel space points
      std::string outputPixelSpacePoints = "outputPixelSpacepoints";
      // name of the output strip space points
      std::string outputStripSpacePoints = "outputStripSpacepoints";
      // name of the output space points
      std::string outputSpacePoints = "output_spacepoints"; 
    };
    
    RootAthInputReader(const RootAthInputReader &) = delete;
    RootAthInputReader(const RootAthInputReader &&) = delete;
    
    
    // Constructor
    /// @param config The configuration struct
    RootAthInputReader(const Config &config, Acts::Logging::Level level);
    
    std::string name() const override {return "RootAthInputReader";}
    
    /// Return the available events range.
    std::pair<std::size_t, std::size_t> availableEvents() const override {
      return {0u,m_events};
    }

    /// Read out data from the input stream
    ///
    /// @param context The algorithm context
    ProcessCode read(const ActsExamples::AlgorithmContext &context) override;

    /// Readonly access to the config
    const Config &config() const { return m_cfg; }
    
  private:
    
    /// Private access to the logging instance
    const Acts::Logger &logger() const {return *m_logger;}
    
    /// The config class
    Config m_cfg;
    
    WriteDataHandle<SimSpacePointContainer> m_outputPixelSpacePoints{this,"outputPixelSpacepoints"};
    WriteDataHandle<SimSpacePointContainer> m_outputStripSpacePoints{this,"outputStripSpacepoints"};
    WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{this, "output_spacepoints"};
    std::unique_ptr<const Acts::Logger> m_logger;
    
    std::mutex m_read_mutex;
    
    /// Vector of {eventNr, entryMin, entryMax}
    std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;
    std::shared_ptr<TChain> m_inputchain;
    long unsigned int m_events;
    
    
    static const unsigned int maxCL  = 1500000;
    static const unsigned int maxSP  = 1500000;
    static const unsigned int maxDTT = 1500000;
    static const unsigned int maxTRK = 1500000;
    static const unsigned int maxPart = 1500000;
    
    // Declaration of leaf types
    unsigned int         run_number;
    long unsigned int    event_number;
    int                  nSE;
    int                  SEID[4];   //[nSE]
    int                  nCL;
    int                  CLindex[maxCL];   //[nCL]
    std::vector<std::string>       *CLhardware;
    double               CLx[maxCL];   //[nCL]
    Double_t        CLy[maxCL];   //[nCL]
    Double_t        CLz[maxCL];   //[nCL]
    Int_t           CLbarrel_endcap[maxCL];   //[nCL]
    Int_t           CLlayer_disk[maxCL];   //[nCL]
    Int_t           CLeta_module[maxCL];   //[nCL]
    Int_t           CLphi_module[maxCL];   //[nCL]
    Int_t           CLside[maxCL];   //[nCL]
    ULong64_t       CLmoduleID[maxCL];   //[nCL]
    std::vector<std::vector<int> > *CLparticleLink_eventIndex;
    std::vector<std::vector<int> > *CLparticleLink_barcode;
    std::vector<std::vector<bool> > *CLbarcodesLinked;
    std::vector<std::vector<float> > *CLparticle_charge;
    std::vector<std::vector<int> > *CLphis;
    std::vector<std::vector<int> > *CLetas;
    std::vector<std::vector<int> > *CLtots;
    Double_t        CLloc_direction1[maxCL];   //[nCL]
    Double_t        CLloc_direction2[maxCL];   //[nCL]
    Double_t        CLloc_direction3[maxCL];   //[nCL]
    Double_t        CLJan_loc_direction1[maxCL];   //[nCL]
    Double_t        CLJan_loc_direction2[maxCL];   //[nCL]
    Double_t        CLJan_loc_direction3[maxCL];   //[nCL]
    Int_t           CLpixel_count[maxCL];   //[nCL]
    Float_t         CLcharge_count[maxCL];   //[nCL]
    Float_t         CLloc_eta[maxCL];   //[nCL]
    Float_t         CLloc_phi[maxCL];   //[nCL]
    Float_t         CLglob_eta[maxCL];   //[nCL]
    Float_t         CLglob_phi[maxCL];   //[nCL]
    Double_t        CLeta_angle[maxCL];   //[nCL]
    Double_t        CLphi_angle[maxCL];   //[nCL]
    Float_t         CLnorm_x[maxCL];   //[nCL]
    Float_t         CLnorm_y[maxCL];   //[nCL]
    Float_t         CLnorm_z[maxCL];   //[nCL]
    std::vector<std::vector<double> > *CLlocal_cov;
    Int_t           nPartEVT;
    Int_t           Part_event_number[maxPart];   //[nPartEVT]
    Int_t           Part_barcode[maxPart];   //[nPartEVT]
    Float_t         Part_px[maxPart];   //[nPartEVT]
    Float_t         Part_py[maxPart];   //[nPartEVT]
    Float_t         Part_pz[maxPart];   //[nPartEVT]
    Float_t         Part_pt[maxPart];   //[nPartEVT]
    Float_t         Part_eta[maxPart];   //[nPartEVT]
    Float_t         Part_vx[maxPart];   //[nPartEVT]
    Float_t         Part_vy[maxPart];   //[nPartEVT]
    Float_t         Part_vz[maxPart];   //[nPartEVT]
    Float_t         Part_radius[maxPart];   //[nPartEVT]
    Float_t         Part_status[maxPart];   //[nPartEVT]
    Float_t         Part_charge[maxPart];   //[nPartEVT]
    Int_t           Part_pdg_id[maxPart];   //[nPartEVT]
    Int_t           Part_passed[maxPart];   //[nPartEVT]
    Int_t           Part_vProdNin[maxPart];   //[nPartEVT]
    Int_t           Part_vProdNout[maxPart];   //[nPartEVT]
    Int_t           Part_vProdStatus[maxPart];   //[nPartEVT]
    Int_t           Part_vProdBarcode[maxPart];   //[nPartEVT]
    std::vector<std::vector<int> > *Part_vParentID;
    std::vector<std::vector<int> > *Part_vParentBarcode;
    Int_t           nSP;
    Int_t           SPindex[maxSP];   //[nSP]
    Double_t        SPx[maxSP];   //[nSP]
    Double_t        SPy[maxSP];   //[nSP]
    Double_t        SPz[maxSP];   //[nSP]
    Int_t           SPCL1_index[maxSP];   //[nSP]
    Int_t           SPCL2_index[maxSP];   //[nSP]
    Int_t           SPisOverlap[maxSP];   //[nSP]
    double          SPradius[maxSP]; //[nSP]
    double          SPcovr[maxSP]; //[nSP]
    double          SPcovz[maxSP]; //[nSP]
    float           SPhl_topstrip[maxSP]; //[nSP]
    float           SPhl_botstrip[maxSP]; //[nSP]
    std::vector<std::vector<float>> *SPtopStripDirection;
    std::vector<std::vector<float>> *SPbottomStripDirection;
    std::vector<std::vector<float>> *SPstripCenterDistance;
    std::vector<std::vector<float>> *SPtopStripCenterPosition;
    Int_t           nTRK;
    Int_t           TRKindex[maxTRK];   //[nTRK]
    Int_t           TRKtrack_fitter[maxTRK];   //[nTRK]
    Int_t           TRKparticle_hypothesis[maxTRK];   //[nTRK]
    std::vector<std::vector<int> > *TRKproperties;
    std::vector<std::vector<int> > *TRKpattern;
    Int_t           TRKndof[maxTRK];   //[nTRK]
    Int_t           TRKmot[maxTRK];   //[nTRK]
    Int_t           TRKoot[maxTRK];   //[nTRK]
    Float_t         TRKchiSq[maxTRK];   //[nTRK]
    std::vector<std::vector<int> > *TRKmeasurementsOnTrack_pixcl_sctcl_index;
    std::vector<std::vector<int> > *TRKoutliersOnTrack_pixcl_sctcl_index;
    Int_t           TRKcharge[maxTRK];   //[nTRK]
    std::vector<std::vector<double> > *TRKperigee_position;
    std::vector<std::vector<double> > *TRKperigee_momentum;
    Int_t           TTCindex[maxTRK];   //[nTRK]
    Int_t           TTCevent_index[maxTRK];   //[nTRK]
    Int_t           TTCparticle_link[maxTRK];   //[nTRK]
    Float_t         TTCprobability[maxTRK];   //[nTRK]
    Int_t           nDTT;
    Int_t           DTTindex[maxDTT];   //[nDTT]
    Int_t           DTTsize[maxDTT];   //[nDTT]
    std::vector<std::vector<int> > *DTTtrajectory_eventindex;
    std::vector<std::vector<int> > *DTTtrajectory_barcode;
    std::vector<std::vector<int> > *DTTstTruth_subDetType;
    std::vector<std::vector<int> > *DTTstTrack_subDetType;
    std::vector<std::vector<int> > *DTTstCommon_subDetType;
    
    
    // List of branches
    TBranch        *b_run_number;   //!
    TBranch        *b_event_number;   //!
    TBranch        *b_nSE;   //!
    TBranch        *b_SEID;   //!
    TBranch        *b_nCL;   //!
    TBranch        *b_CLindex;   //!
    TBranch        *b_CLhardware;   //!
    TBranch        *b_CLx;   //!
    TBranch        *b_CLy;   //!
    TBranch        *b_CLz;   //!
    TBranch        *b_CLbarrel_endcap;   //!
    TBranch        *b_CLlayer_disk;   //!
    TBranch        *b_CLeta_module;   //!
    TBranch        *b_CLphi_module;   //!
    TBranch        *b_CLside;   //!
    TBranch        *b_CLmoduleID;   //!
    TBranch        *b_CLparticleLink_eventIndex;   //!
    TBranch        *b_CLparticleLink_barcode;   //!
    TBranch        *b_CLbarcodesLinked;   //!
    TBranch        *b_CLparticle_charge;   //!
    TBranch        *b_CLphis;   //!
    TBranch        *b_CLetas;   //!
    TBranch        *b_CLtots;   //!
    TBranch        *b_CLloc_direction1;   //!
    TBranch        *b_CLloc_direction2;   //!
    TBranch        *b_CLloc_direction3;   //!
    TBranch        *b_CLJan_loc_direction1;   //!
    TBranch        *b_CLJan_loc_direction2;   //!
    TBranch        *b_CLJan_loc_direction3;   //!
    TBranch        *b_CLpixel_count;   //!
    TBranch        *b_CLcharge_count;   //!
    TBranch        *b_CLloc_eta;   //!
    TBranch        *b_CLloc_phi;   //!
    TBranch        *b_CLglob_eta;   //!
    TBranch        *b_CLglob_phi;   //!
    TBranch        *b_CLeta_angle;   //!
    TBranch        *b_CLphi_angle;   //!
    TBranch        *b_CLnorm_x;   //!
    TBranch        *b_CLnorm_y;   //!
    TBranch        *b_CLnorm_z;   //!
    TBranch        *b_CLlocal_cov;   //!
    TBranch        *b_nPartEVT;   //!
    TBranch        *b_Part_event_number;   //!
    TBranch        *b_Part_barcode;   //!
    TBranch        *b_Part_px;   //!
    TBranch        *b_Part_py;   //!
    TBranch        *b_Part_pz;   //!
    TBranch        *b_Part_pt;   //!
    TBranch        *b_Part_eta;   //!
    TBranch        *b_Part_vx;   //!
    TBranch        *b_Part_vy;   //!
    TBranch        *b_Part_vz;   //!
    TBranch        *b_Part_radius;   //!
    TBranch        *b_Part_status;   //!
    TBranch        *b_Part_charge;   //!
    TBranch        *b_Part_pdg_id;   //!
    TBranch        *b_Part_passed;   //!
    TBranch        *b_Part_vProdNin;   //!
    TBranch        *b_Part_vProdNout;   //!
    TBranch        *b_Part_vProdStatus;   //!
    TBranch        *b_Part_vProdBarcode;   //!
    TBranch        *b_Part_vParentID;   //!
    TBranch        *b_Part_vParentBarcode;   //!
    TBranch        *b_nSP;   //!
    TBranch        *b_SPindex;   //!
    TBranch        *b_SPx;   //!
    TBranch        *b_SPy;   //!
    TBranch        *b_SPz;   //!
    TBranch        *b_SPCL1_index;   //!
    TBranch        *b_SPCL2_index;   //!
    TBranch        *b_SPisOverlap;   //!
    TBranch        *b_SPradius; //!
    TBranch        *b_SPcovr; //!
    TBranch        *b_SPcovz; //!
    TBranch        *b_SPhl_topstrip; //!
    TBranch        *b_SPhl_botstrip; //!
    TBranch        *b_SPtopStripDirection; //!
    TBranch        *b_SPbottomStripDirection; //!
    TBranch        *b_SPstripCenterDistance; //!
    TBranch        *b_SPtopStripCenterPosition;//!
    TBranch        *b_nTRK;   //!
    TBranch        *b_TRKindex;   //!
    TBranch        *b_TRKtrack_fitter;   //!
    TBranch        *b_TRKparticle_hypothesis;   //!
    TBranch        *b_TRKproperties;   //!
    TBranch        *b_TRKpattern;   //!
    TBranch        *b_TRKndof;   //!
    TBranch        *b_TRKmot;   //!
    TBranch        *b_TRKoot;   //!
    TBranch        *b_TRKchiSq;   //!
    TBranch        *b_TRKmeasurementsOnTrack_pixcl_sctcl_index;   //!
    TBranch        *b_TRKoutliersOnTrack_pixcl_sctcl_index;   //!
    TBranch        *b_TRKcharge;   //!
    TBranch        *b_TRKperigee_position;   //!
    TBranch        *b_TRKperigee_momentum;   //!
    TBranch        *b_TTCindex;   //!
    TBranch        *b_TTCevent_index;   //!
    TBranch        *b_TTCparticle_link;   //!
    TBranch        *b_TTCprobability;   //!
    TBranch        *b_nDTT;   //!
    TBranch        *b_DTTindex;   //!
    TBranch        *b_DTTsize;   //!
    TBranch        *b_DTTtrajectory_eventindex;   //!
    TBranch        *b_DTTtrajectory_barcode;   //!
    TBranch        *b_DTTstTruth_subDetType;   //!
    TBranch        *b_DTTstTrack_subDetType;   //!
    TBranch        *b_DTTstCommon_subDetType;   //!
    
    
    
    
    
    
    
    
  };
} // ActsExamples
