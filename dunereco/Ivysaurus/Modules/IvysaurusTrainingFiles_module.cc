/**
 *  @file   dunereco/Ivysaurus/Modules/IvysaurusTrainingFiles_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "dunereco/Ivysaurus/Managers/GridManager.h"
#include "dunereco/Ivysaurus/Managers/TrackVarManager.h"
#include "dunereco/Ivysaurus/Managers/ShowerVarManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"
#include "dunereco/Ivysaurus/TensorFlow/IvysaurusEvaluator.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ivysaurus
{

/**
 *  @brief  IvysaurusTrainingFiles class
 */
class IvysaurusTrainingFiles : public art::EDAnalyzer
{
public:
   IvysaurusTrainingFiles(fhicl::ParameterSet const &pset);
   virtual ~IvysaurusTrainingFiles();

   void beginJob();
   void endJob();
   void analyze(const art::Event &evt);

   void Reset();

private:

  // Trees
  TTree *m_tree;

  // Tree variables
  int m_run;
  int m_subrun;
  int m_event;

  int m_truePDG;
  float m_completeness;
  float m_purity;

  int m_nSpacePoints;
  int m_nHits2D;

  float m_trackScore;

  std::vector<std::vector<double>> m_spacePoints;
  std::vector<std::vector<double>> m_projectionsU;
  std::vector<std::vector<double>> m_projectionsV;
  std::vector<std::vector<double>> m_projectionsW;
  std::vector<float> m_startDriftBoundariesU;
  std::vector<float> m_startDriftBoundariesV;
  std::vector<float> m_startDriftBoundariesW;
  std::vector<float> m_endDriftBoundariesU;
  std::vector<float> m_endDriftBoundariesV;
  std::vector<float> m_endDriftBoundariesW;
  std::vector<float> m_startWireBoundariesU;
  std::vector<float> m_startWireBoundariesV;
  std::vector<float> m_startWireBoundariesW;
  std::vector<float> m_endWireBoundariesU;
  std::vector<float> m_endWireBoundariesV;
  std::vector<float> m_endWireBoundariesW;
  std::vector<std::vector<float>> m_startGridValuesU;
  std::vector<std::vector<float>> m_startGridValuesV;
  std::vector<std::vector<float>> m_startGridValuesW;
  std::vector<std::vector<float>> m_endGridValuesU;
  std::vector<std::vector<float>> m_endGridValuesV;
  std::vector<std::vector<float>> m_endGridValuesW;
  std::vector<std::vector<float>> m_startGridValuesU_disp;
  std::vector<std::vector<float>> m_startGridValuesV_disp;
  std::vector<std::vector<float>> m_startGridValuesW_disp;
  std::vector<std::vector<float>> m_endGridValuesU_disp;
  std::vector<std::vector<float>> m_endGridValuesV_disp;
  std::vector<std::vector<float>> m_endGridValuesW_disp;

  // TrackVars
  int m_trackVarsSuccessful;
  float m_nTrackChildren;
  float m_nShowerChildren;
  float m_nGrandChildren;
  float m_nChildHits;
  float m_childEnergy;
  float m_childTrackScore;
  float m_trackLength;
  float m_trackWobble;
  float m_trackMomComparison;

  // ShowerVars
  int m_showerVarsSuccessful;
  float m_showerDisplacement;
  float m_DCA;
  float m_trackStubLength;
  float m_nuVertexAvSeparation;
  float m_nuVertexChargeAsymmetry;
  float m_showerFoundConnectionPathway;
  float m_showerInitialGapSize;
  float m_showerLargestGapSize;
  float m_showerPathwayLength;
  float m_showerPathwayScatteringAngle2D;
  float m_showerNHits;
  float m_showerFoundHitRatio;
  float m_showerScatterAngle;
  float m_showerOpeningAngle;
  float m_showerNuVertexEnergyAsymmetry;
  float m_showerNuVertexEnergyWeightedMeanRadialDistance;
  float m_showerStartEnergyAsymmetry;
  float m_showerStartMoliereRadius;
  float m_showerNAmbiguousViews;
  float m_showerUnaccountedEnergy;

  // Managers
  GridManager m_gridManager;
  TrackVarManager m_trackVarManager;
  ShowerVarManager m_showerVarManager;
  //IvysaurusEvaluator m_ivysaurusEvaluator;

  // FCL module labels
  std::string m_hitModuleLabel;
  std::string m_recoModuleLabel;

  // Module variables
  float m_completenessThreshold;
  float m_purityThreshold;
  float m_nSpacepointThreshold;
  bool m_writeVisualisationInfo;
};

DEFINE_ART_MODULE(IvysaurusTrainingFiles)

} // namespace ivysaurus

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include <iostream>
#include <random>

namespace ivysaurus
{

IvysaurusTrainingFiles::IvysaurusTrainingFiles(fhicl::ParameterSet const &pset) : 
    art::EDAnalyzer(pset),
    m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_showerVarManager(pset.get<fhicl::ParameterSet>("ShowerVarManager")),
    //m_ivysaurusEvaluator(pset.get<fhicl::ParameterSet>("IvysaurusEvaluator")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_completenessThreshold(pset.get<float>("CompletenessThreshold")),
    m_purityThreshold(pset.get<float>("PurityThreshold")),
    m_nSpacepointThreshold(pset.get<float>("NSpacepointThreshold")),
    m_writeVisualisationInfo(pset.get<bool>("WriteVisualisationInfo"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

IvysaurusTrainingFiles::~IvysaurusTrainingFiles()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::analyze(const art::Event &evt)
{
    std::cout << "BEGIN" << std::endl;

    const std::vector<art::Ptr<recob::PFParticle>> pfparticles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, m_recoModuleLabel);

    // Get the neutrino PFP
    art::Ptr<recob::PFParticle> nuPFP;

    try
    {
        nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    }
    catch(...)
    {
        return;
    }

    const std::vector<art::Ptr<recob::PFParticle>> &nuChildPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, m_recoModuleLabel);

    for (const art::Ptr<recob::PFParticle> &pfparticle : pfparticles)
    {
        bool isPrimary = false;

        for (const art::Ptr<recob::PFParticle> &nuChildPFP : nuChildPFPs)
        {
            if (nuChildPFP->Self() == pfparticle->Self())
            {
                isPrimary = true;
                break;
            }
        }

        if (!isPrimary)
            continue;

        Reset();

        m_run = evt.run();
        m_subrun = evt.subRun();
        m_event = evt.event();

        ////////////////////////////////////////////                                                                                                                                                                                    
        // First, let's get the truth information...
        ////////////////////////////////////////////  
        const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);
        const std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, m_hitModuleLabel);

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);

        if (TruthMatchUtils::Valid(g4id))
        {
            art::ServiceHandle<cheat::ParticleInventoryService> piServ;
            m_completeness = IvysaurusUtils::CompletenessFromTrueParticleID(clockData, pfpHits, eventHitList, g4id);
            m_purity = IvysaurusUtils::HitPurityFromTrueParticleID(clockData, pfpHits, g4id);
            m_truePDG = piServ->ParticleList().at(g4id)->PdgCode();
        }
        else
        {
            continue;
        }

        // If it isn't a PDG that we care about, move on...
        int absPDG = std::abs(m_truePDG);

        if ((absPDG != 13) && (absPDG != 2212) && (absPDG != 211) && (absPDG != 11) & (absPDG != 22))
            continue;

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now, get the track score... 
        ////////////////////////////////////////////  
        const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
        const auto metaMap = metadata->GetPropertiesMap();

        if (metaMap.find("TrackScore") == metaMap.end())
            continue;

        m_trackScore = metaMap.at("TrackScore");

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now, get space points into file.. 
        ////////////////////////////////////////////  
        const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);

        if (spacepoints.empty())
            continue;

        for (art::Ptr<recob::SpacePoint> spacepoint : spacepoints)
            m_spacePoints.push_back({spacepoint->XYZ()[0], spacepoint->XYZ()[1], spacepoint->XYZ()[2]});

        m_nSpacePoints = spacepoints.size();
        m_nHits2D = pfpHits.size();

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now, apply quality cuts
        ////////////////////////////////////////////  

        //std::cout << "m_completeness: " << m_completeness << std::endl;
        //std::cout << "m_purity: " << m_purity << std::endl;
        //std::cout << "spacepoints.size(): " << spacepoints << std::endl;

        if ((m_completeness < m_completenessThreshold) || (m_purity < m_purityThreshold) || (spacepoints.size() < m_nSpacepointThreshold))
            continue;

        ////////////////////////////////////////////
        // Now, into 2D
        ////////////////////////////////////////////  

        int nInitialisedGrids(0);

        for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U, IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
        {
            GridManager::Grid startGrid = m_gridManager.ObtainViewGrid(evt, pfparticle, pandoraView, true);
            GridManager::Grid endGrid = m_gridManager.ObtainViewGrid(evt, pfparticle, pandoraView, false);

            if (!startGrid.IsInitialised() || !endGrid.IsInitialised())
                continue;

            ++nInitialisedGrids;

            ////////////////////////////////////////////
            // Let's fill the hit vectors...
            ////////////////////////////////////////////
            std::vector<std::vector<double>> &projections = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_projectionsU :
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_projectionsV : m_projectionsW;

            for (const art::Ptr<recob::Hit> hit : pfpHits)
            {
                // Make sure 2D hit has an associated space point
                std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

                if (spacePoints.empty())
                    continue;

                const IvysaurusUtils::PandoraView thisPandoraView = IvysaurusUtils::GetPandora2DView(hit);

                if (thisPandoraView != pandoraView)
                    continue;

                // Get 2D hit position
                float width = 0.f;
                TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
                IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, pandoraView, pandoraHitPosition, width);
                
                if (m_writeVisualisationInfo)
                {
                    projections.push_back({pandoraHitPosition.X(), pandoraHitPosition.Y(), pandoraHitPosition.Z()});
                }
                else
                {
                    if (startGrid.IsInsideGrid(pandoraHitPosition, width) || endGrid.IsInsideGrid(pandoraHitPosition, width))
                        projections.push_back({pandoraHitPosition.X(), pandoraHitPosition.Y(), pandoraHitPosition.Z()});
                }
            }

            ////////////////////////////////////////////
            // Then fill grid vectors...
            ////////////////////////////////////////////
            // CaloGrid
            m_gridManager.FillViewGrid(evt, pfparticle, startGrid);
            m_gridManager.FillViewGrid(evt, pfparticle, endGrid);

            // DisplacementGrid
            const art::Ptr<recob::Vertex> nuVertex3D = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
            const TVector3 nuVertex3D_tv = TVector3(nuVertex3D->position().X(), nuVertex3D->position().Y(), nuVertex3D->position().Z());
            GridManager::Grid startGrid_disp = m_gridManager.ObtainViewDisplacementGrid(evt, nuVertex3D_tv, startGrid);
            GridManager::Grid endGrid_disp = m_gridManager.ObtainViewDisplacementGrid(evt, nuVertex3D_tv, endGrid);

            std::vector<float> &startDriftBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startDriftBoundariesU :
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startDriftBoundariesV : m_startDriftBoundariesW;

            std::vector<float> &endDriftBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endDriftBoundariesU :
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endDriftBoundariesV : m_endDriftBoundariesW;

            std::vector<float> &startWireBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startWireBoundariesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startWireBoundariesV : m_startWireBoundariesW;

            std::vector<float> &endWireBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endWireBoundariesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endWireBoundariesV : m_endWireBoundariesW;

            std::vector<std::vector<float>> &startGridValues = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startGridValuesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startGridValuesV : m_startGridValuesW;

            std::vector<std::vector<float>> &endGridValues = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endGridValuesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endGridValuesV : m_endGridValuesW;

            std::vector<std::vector<float>> &startGridValues_disp = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startGridValuesU_disp : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startGridValuesV_disp : m_startGridValuesW_disp;

            std::vector<std::vector<float>> &endGridValues_disp = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endGridValuesU_disp : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endGridValuesV_disp : m_endGridValuesW_disp;

            startDriftBoundaries = startGrid.GetDriftBoundaries();
            endDriftBoundaries = endGrid.GetDriftBoundaries();
            startWireBoundaries = startGrid.GetWireBoundaries();
            endWireBoundaries = endGrid.GetWireBoundaries();
            startGridValues = startGrid.GetGridValues();
            endGridValues = endGrid.GetGridValues();
            startGridValues_disp = startGrid_disp.GetGridValues();
            endGridValues_disp = endGrid_disp.GetGridValues();
        }

        if (nInitialisedGrids != 3)
            continue;

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Just leave this here a minute
        ////////////////////////////////////////////  
        //m_ivysaurusEvaluator.IvysaurusUseEvaluate(evt, pfparticle);

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now fill the track variables
        ////////////////////////////////////////////  
        TrackVarManager::TrackVars trackVars;

        m_trackVarsSuccessful = m_trackVarManager.EvaluateTrackVars(evt, pfparticle, trackVars) ? 1 : 0;
        m_trackVarManager.NormaliseTrackVars(trackVars);

        m_nTrackChildren = trackVars.GetNTrackChildren();
        m_nShowerChildren = trackVars.GetNShowerChildren();
        m_nGrandChildren = trackVars.GetNGrandChildren();
        m_nChildHits = trackVars.GetNChildHits();
        m_childEnergy = trackVars.GetChildEnergy();
        m_childTrackScore = trackVars.GetChildTrackScore();
        m_trackLength = trackVars.GetTrackLength();
        m_trackWobble = trackVars.GetWobble();
        m_trackMomComparison = trackVars.GetMomentumComparison();

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now fill the shower variables
        ////////////////////////////////////////////  
        ShowerVarManager::ShowerVars showerVars;

        m_showerVarsSuccessful = m_showerVarManager.EvaluateShowerVars(evt, pfparticle, showerVars) ? 1 : 0;
        m_showerVarManager.NormaliseShowerVars(showerVars);

        m_showerDisplacement = showerVars.GetDisplacement();
        m_DCA = showerVars.GetDCA();
        m_trackStubLength = showerVars.GetTrackStubLength();
        m_nuVertexAvSeparation = showerVars.GetNuVertexAvSeparation();
        m_nuVertexChargeAsymmetry = showerVars.GetNuVertexChargeAsymmetry();
        m_showerFoundConnectionPathway = showerVars.GetFoundConnectionPathway();
        m_showerInitialGapSize = showerVars.GetInitialGapSize();
        m_showerLargestGapSize = showerVars.GetLargestGapSize();
        m_showerPathwayLength = showerVars.GetPathwayLength();
        m_showerPathwayScatteringAngle2D = showerVars.GetPathwayScatteringAngle2D();
        m_showerNHits = showerVars.GetNShowerHits();
        m_showerFoundHitRatio = showerVars.GetFoundHitRatio();
        m_showerScatterAngle = showerVars.GetScatterAngle();
        m_showerOpeningAngle = showerVars.GetOpeningAngle();
        m_showerNuVertexEnergyAsymmetry = showerVars.GetNuVertexEnergyAsymmetry();
        m_showerNuVertexEnergyWeightedMeanRadialDistance = showerVars.GetNuVertexEnergyWeightedMeanRadialDistance();
        m_showerStartEnergyAsymmetry = showerVars.GetShowerStartEnergyAsymmetry();
        m_showerStartMoliereRadius = showerVars.GetShowerStartMoliereRadius();
        m_showerNAmbiguousViews = showerVars.GetNAmbiguousViews();
        m_showerUnaccountedEnergy = showerVars.GetUnaccountedEnergy();

        m_tree->Fill();
    }
    std::cout << "DONE" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::Reset()
{
  const int defaultInt = -999;
  const float defaultFloat = -999.f;

  m_run = defaultInt;
  m_subrun = defaultInt;
  m_event = defaultInt;

  m_truePDG = defaultInt;
  m_completeness = defaultFloat;
  m_purity = defaultFloat;

  m_nSpacePoints = defaultInt;
  m_nHits2D = defaultInt;
  m_trackScore = defaultFloat;

  m_spacePoints.clear();
  m_projectionsU.clear();
  m_projectionsV.clear();
  m_projectionsW.clear();
  m_startDriftBoundariesU.clear();
  m_startDriftBoundariesV.clear();
  m_startDriftBoundariesW.clear();
  m_endDriftBoundariesU.clear();
  m_endDriftBoundariesV.clear();
  m_endDriftBoundariesW.clear();
  m_startWireBoundariesU.clear();
  m_startWireBoundariesV.clear();
  m_startWireBoundariesW.clear();
  m_endWireBoundariesU.clear();
  m_endWireBoundariesV.clear();
  m_endWireBoundariesW.clear();
  m_startGridValuesU.clear();
  m_startGridValuesV.clear();
  m_startGridValuesW.clear();
  m_endGridValuesU.clear();
  m_endGridValuesV.clear();
  m_endGridValuesW.clear();
  m_startGridValuesU_disp.clear();
  m_startGridValuesV_disp.clear();
  m_startGridValuesW_disp.clear();
  m_endGridValuesU_disp.clear();
  m_endGridValuesV_disp.clear();
  m_endGridValuesW_disp.clear();

  m_trackVarsSuccessful = 0;
  m_nTrackChildren = defaultFloat;
  m_nShowerChildren = defaultFloat;
  m_nGrandChildren = defaultFloat;
  m_nChildHits = defaultFloat;
  m_childEnergy = defaultFloat;
  m_childTrackScore = defaultFloat;
  m_trackLength = defaultFloat;
  m_trackWobble = defaultFloat;
  m_trackMomComparison = defaultFloat;

  m_showerVarsSuccessful = 0;
  m_showerDisplacement = defaultFloat;
  m_DCA = defaultFloat;
  m_trackStubLength = defaultFloat;
  m_nuVertexAvSeparation = defaultFloat;
  m_nuVertexChargeAsymmetry = defaultFloat;
  m_showerFoundConnectionPathway = defaultFloat;
  m_showerInitialGapSize = defaultFloat;
  m_showerLargestGapSize = defaultFloat;
  m_showerPathwayLength = defaultFloat;
  m_showerPathwayScatteringAngle2D = defaultFloat;
  m_showerNHits = defaultFloat;
  m_showerFoundHitRatio = defaultFloat;
  m_showerScatterAngle = defaultFloat;
  m_showerOpeningAngle = defaultFloat;
  m_showerNuVertexEnergyAsymmetry = defaultFloat;
  m_showerNuVertexEnergyWeightedMeanRadialDistance = defaultFloat;
  m_showerStartEnergyAsymmetry = defaultFloat;
  m_showerStartMoliereRadius = defaultFloat;
  m_showerNAmbiguousViews = defaultFloat;
  m_showerUnaccountedEnergy = defaultFloat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::beginJob()
{
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    m_tree = tfs->make<TTree>("ivysaur", "Ivysaur");

    m_tree->Branch("Run", &m_run);
    m_tree->Branch("Subrun", &m_subrun);
    m_tree->Branch("Event", &m_event);

    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("Completeness", &m_completeness);
    m_tree->Branch("Purity", &m_purity);

    m_tree->Branch("NSpacePoints", &m_nSpacePoints);
    m_tree->Branch("NHits2D", &m_nHits2D);

    m_tree->Branch("TrackScore", &m_trackScore);

    m_tree->Branch("SpacePoints", &m_spacePoints);
    m_tree->Branch("ProjectionsU", &m_projectionsU);
    m_tree->Branch("ProjectionsV", &m_projectionsV);
    m_tree->Branch("ProjectionsW", &m_projectionsW);

    m_tree->Branch("StartDriftBoundariesU", &m_startDriftBoundariesU);
    m_tree->Branch("StartDriftBoundariesV", &m_startDriftBoundariesV);
    m_tree->Branch("StartDriftBoundariesW", &m_startDriftBoundariesW);
    m_tree->Branch("EndDriftBoundariesU", &m_endDriftBoundariesU);
    m_tree->Branch("EndDriftBoundariesV", &m_endDriftBoundariesV);
    m_tree->Branch("EndDriftBoundariesW", &m_endDriftBoundariesW);
    m_tree->Branch("StartWireBoundariesU", &m_startWireBoundariesU);
    m_tree->Branch("StartWireBoundariesV", &m_startWireBoundariesV);
    m_tree->Branch("StartWireBoundariesW", &m_startWireBoundariesW);
    m_tree->Branch("EndWireBoundariesU", &m_endWireBoundariesU);
    m_tree->Branch("EndWireBoundariesV", &m_endWireBoundariesV);
    m_tree->Branch("EndWireBoundariesW", &m_endWireBoundariesW);
    m_tree->Branch("StartGridU", &m_startGridValuesU);
    m_tree->Branch("StartGridV", &m_startGridValuesV);
    m_tree->Branch("StartGridW", &m_startGridValuesW);
    m_tree->Branch("EndGridU", &m_endGridValuesU);
    m_tree->Branch("EndGridV", &m_endGridValuesV);
    m_tree->Branch("EndGridW", &m_endGridValuesW);
    m_tree->Branch("StartGridUDisp", &m_startGridValuesU_disp);
    m_tree->Branch("StartGridVDisp", &m_startGridValuesV_disp);
    m_tree->Branch("StartGridWDisp", &m_startGridValuesW_disp);
    m_tree->Branch("EndGridUDisp", &m_endGridValuesU_disp);
    m_tree->Branch("EndGridVDisp", &m_endGridValuesV_disp);
    m_tree->Branch("EndGridWDisp", &m_endGridValuesW_disp);

    m_tree->Branch("TrackVarsSuccessful", &m_trackVarsSuccessful);
    m_tree->Branch("NTrackChildren", &m_nTrackChildren);
    m_tree->Branch("NShowerChildren", &m_nShowerChildren);
    m_tree->Branch("NGrandChildren", &m_nGrandChildren);
    m_tree->Branch("NChildHits", &m_nChildHits);
    m_tree->Branch("ChildEnergy", &m_childEnergy);
    m_tree->Branch("ChildTrackScore", &m_childTrackScore);
    m_tree->Branch("TrackLength", &m_trackLength);
    m_tree->Branch("TrackWobble", &m_trackWobble);
    m_tree->Branch("TrackMomComparison", &m_trackMomComparison);

    m_tree->Branch("ShowerVarsSuccessful", &m_showerVarsSuccessful);
    m_tree->Branch("ShowerDisplacement", &m_showerDisplacement);
    m_tree->Branch("ShowerDCA", &m_DCA);
    m_tree->Branch("ShowerTrackStubLength", &m_trackStubLength);
    m_tree->Branch("ShowerNuVertexAvSeparation", &m_nuVertexAvSeparation);
    m_tree->Branch("ShowerNuVertexChargeAsymmetry", &m_nuVertexChargeAsymmetry);
    m_tree->Branch("ShowerFoundConnectionPathway", &m_showerFoundConnectionPathway);
    m_tree->Branch("ShowerInitialGapSize", &m_showerInitialGapSize);
    m_tree->Branch("ShowerLargestGapSize", &m_showerLargestGapSize);
    m_tree->Branch("ShowerPathwayLength", &m_showerPathwayLength);
    m_tree->Branch("ShowerPathwayScatteringAngle2D", &m_showerPathwayScatteringAngle2D);
    m_tree->Branch("ShowerNHits", &m_showerNHits);
    m_tree->Branch("ShowerFoundHitRatio", &m_showerFoundHitRatio);
    m_tree->Branch("ShowerScatterAngle", &m_showerScatterAngle);
    m_tree->Branch("ShowerOpeningAngle", &m_showerOpeningAngle);
    m_tree->Branch("ShowerNuVertexEnergyAsymmetry", &m_showerNuVertexEnergyAsymmetry);
    m_tree->Branch("ShowerNuVertexEnergyWeightedMeanRadialDistance", &m_showerNuVertexEnergyWeightedMeanRadialDistance);
    m_tree->Branch("ShowerStartEnergyAsymmetry", &m_showerStartEnergyAsymmetry);
    m_tree->Branch("ShowerStartMoliereRadius", &m_showerStartMoliereRadius);
    m_tree->Branch("ShowerNAmbiguousViews", &m_showerNAmbiguousViews);
    m_tree->Branch("ShowerUnaccountedEnergy", &m_showerUnaccountedEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::endJob()
{
}


} //namespace ivysaurus

