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

  std::vector<int> m_truePDG;
  std::vector<float> m_completeness;
  std::vector<float> m_purity;

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

  // Managers
  GridManager m_gridManager;

  // FCL module labels
  std::string m_hitModuleLabel;
  std::string m_recoModuleLabel;
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
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

IvysaurusTrainingFiles::~IvysaurusTrainingFiles()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::analyze(const art::Event &evt)
{
    // Now our grid!
    const std::vector<art::Ptr<recob::PFParticle>> pfparticles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, m_recoModuleLabel);

    for (const art::Ptr<recob::PFParticle> pfparticle : pfparticles)
    {
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
            m_completeness.push_back(IvysaurusUtils::CompletenessFromTrueParticleID(clockData, pfpHits, eventHitList, g4id));
            m_purity.push_back(IvysaurusUtils::HitPurityFromTrueParticleID(clockData, pfpHits, g4id));
            m_truePDG.push_back(piServ->ParticleList().at(g4id)->PdgCode());
        }

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now, get space points into file.. 
        ////////////////////////////////////////////  
        const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);

        for (art::Ptr<recob::SpacePoint> spacepoint : spacepoints)
            m_spacePoints.push_back({spacepoint->XYZ()[0], spacepoint->XYZ()[1], spacepoint->XYZ()[2]});

        ////////////////////////////////////////////                                                                                                                                                                                    
        // Now, into 2D
        ////////////////////////////////////////////  
        for (GridManager::PandoraView pandoraView : {GridManager::PandoraView::TPC_VIEW_U, GridManager::PandoraView::TPC_VIEW_V, GridManager::PandoraView::TPC_VIEW_W})
        {
            GridManager::Grid startGrid = m_gridManager.ObtainViewGrid(evt, pfparticle, pandoraView, true);
            GridManager::Grid endGrid = m_gridManager.ObtainViewGrid(evt, pfparticle, pandoraView, false);

            if (!startGrid.IsInitialised() || !endGrid.IsInitialised())
                continue;

            ////////////////////////////////////////////
            // Let's fill the hit vectors...
            ////////////////////////////////////////////
            std::vector<std::vector<double>> &projections = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_projectionsU :
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_projectionsV : m_projectionsW;

            for (const art::Ptr<recob::Hit> hit : pfpHits)
            {
                // Make sure 2D hit has an associated space point
                std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

                if (spacePoints.empty())
                    continue;
                     
                const GridManager::PandoraView thisPandoraView = m_gridManager.GetPandora2DView(hit);

                if (thisPandoraView != pandoraView)
                    continue;

                // Get 2D hit position
                const TVector3 pandoraHitPosition = m_gridManager.ObtainPandoraHitPosition(evt, hit, pandoraView);

                projections.push_back({pandoraHitPosition.X(), pandoraHitPosition.Y(), pandoraHitPosition.Z()});
            }

            ////////////////////////////////////////////
            // Then fill grid vectors...
            ////////////////////////////////////////////
            m_gridManager.FillViewGrid(evt, pfparticle, startGrid);
            m_gridManager.FillViewGrid(evt, pfparticle, endGrid);

            std::vector<float> &startDriftBoundaries = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_startDriftBoundariesU :
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_startDriftBoundariesV : m_startDriftBoundariesW;

            std::vector<float> &endDriftBoundaries = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_endDriftBoundariesU :
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_endDriftBoundariesV : m_endDriftBoundariesW;

            std::vector<float> &startWireBoundaries = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_startWireBoundariesU : 
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_startWireBoundariesV : m_startWireBoundariesW;

            std::vector<float> &endWireBoundaries = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_endWireBoundariesU : 
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_endWireBoundariesV : m_endWireBoundariesW;

            std::vector<std::vector<float>> &startGridValues = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_startGridValuesU : 
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_startGridValuesV : m_startGridValuesW;

            std::vector<std::vector<float>> &endGridValues = pandoraView == GridManager::PandoraView::TPC_VIEW_U ? m_endGridValuesU : 
                pandoraView == GridManager::PandoraView::TPC_VIEW_V ? m_endGridValuesV : m_endGridValuesW;

            startDriftBoundaries = startGrid.GetDriftBoundaries();
            endDriftBoundaries = endGrid.GetDriftBoundaries();
            startWireBoundaries = startGrid.GetWireBoundaries();
            endWireBoundaries = endGrid.GetWireBoundaries();
            startGridValues = startGrid.GetGridValues();
            endGridValues = endGrid.GetGridValues();
        }

        m_tree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::Reset()
{
  const int defaultInt = -999;

  m_run = defaultInt;
  m_subrun = defaultInt;
  m_event = defaultInt;

  m_truePDG.clear();
  m_completeness.clear();
  m_purity.clear();

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::endJob()
{
}


} //namespace ivysaurus

