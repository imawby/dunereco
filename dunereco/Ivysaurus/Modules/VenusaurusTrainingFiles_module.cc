/**
 *  @file   dunereco/Ivysaurus/Modules/VenusaurusTrainingFiles_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "dunereco/Ivysaurus/Managers/ProfileManager.h"
#include "dunereco/Ivysaurus/Managers/TrackVarManager.h"
#include "dunereco/Ivysaurus/Managers/ShowerVarManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ivysaurus
{

/**
 *  @brief  VenusaurusTrainingFiles class
 */
class VenusaurusTrainingFiles : public art::EDAnalyzer
{
public:
   VenusaurusTrainingFiles(fhicl::ParameterSet const &pset);
   virtual ~VenusaurusTrainingFiles();

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

  std::vector<float> m_longitudinal_true;
  std::vector<float> m_transverse_true;
  std::vector<float> m_energy_true; 
  std::vector<float> m_fitDirection_true;
  std::vector<std::vector<float>> m_fitPositions_true;

  std::vector<float> m_longitudinal_recoStart;
  std::vector<float> m_transverse_recoStart;
  std::vector<float> m_energy_recoStart; 
  std::vector<float> m_fitDirection_recoStart;
  std::vector<std::vector<float>> m_fitPositions_recoStart;

  std::vector<float> m_longitudinal_pca;
  std::vector<float> m_transverse_pca;
  std::vector<float> m_energy_pca; 
  std::vector<float> m_fitDirection_pca;
  std::vector<std::vector<float>> m_fitPositions_pca;

  std::vector<float> m_longitudinal_trackStub;
  std::vector<float> m_transverse_trackStub;
  std::vector<float> m_energy_trackStub; 
  std::vector<float> m_fitDirection_trackStub;
  std::vector<std::vector<float>> m_fitPositions_trackStub;

  std::vector<float> m_longitudinal_track;
  std::vector<float> m_transverse_track;
  std::vector<float> m_energy_track; 
  std::vector<std::vector<float>> m_fitPositions_track;

  std::vector<std::vector<float>> m_longitudinalProfile;
  std::vector<std::vector<float>> m_transverseProfile;

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
  ProfileManager m_profileManager;
  TrackVarManager m_trackVarManager;
  ShowerVarManager m_showerVarManager;

  // FCL module labels
  std::string m_hitModuleLabel;
  std::string m_recoModuleLabel;

  // Module variables
  float m_completenessThreshold;
  float m_purityThreshold;
  float m_nSpacepointThreshold;
  bool m_writeVisualisationInfo;
};

DEFINE_ART_MODULE(VenusaurusTrainingFiles)

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

VenusaurusTrainingFiles::VenusaurusTrainingFiles(fhicl::ParameterSet const &pset) : 
    art::EDAnalyzer(pset),
    m_profileManager(pset.get<fhicl::ParameterSet>("ProfileManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_showerVarManager(pset.get<fhicl::ParameterSet>("ShowerVarManager")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_completenessThreshold(pset.get<float>("CompletenessThreshold")),
    m_purityThreshold(pset.get<float>("PurityThreshold")),
    m_nSpacepointThreshold(pset.get<float>("NSpacepointThreshold"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VenusaurusTrainingFiles::~VenusaurusTrainingFiles()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VenusaurusTrainingFiles::analyze(const art::Event &evt)
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
        // Now the transverse profiles
        ////////////////////////////////////////////  
        // True method
        TVector3 truePosition = TVector3(0.f, 0.f, 0.f);
        TVector3 trueDirection = TVector3(0.f, 0.f, 0.f);

        if (m_profileManager.GetTrueDirectionAndPosition(evt, pfparticle, trueDirection, truePosition))
        {
            ProfileManager::EnergyProfiles energyProfiles_true = m_profileManager.CreateProfiles(evt, pfparticle, trueDirection, truePosition);
            m_longitudinal_true = energyProfiles_true.m_longitudinal;
            m_transverse_true = energyProfiles_true.m_transverse;
            m_energy_true = energyProfiles_true.m_energy;
            m_fitDirection_true = std::vector<float>({static_cast<float>(trueDirection.X()), static_cast<float>(trueDirection.Y()), static_cast<float>(trueDirection.Z())});
            m_profileManager.FillFitPositions(evt, pfparticle, trueDirection, truePosition, m_fitPositions_true);
        }

        // RecoStart method
        TVector3 recoPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 recoStartDirection = TVector3(0.f, 0.f, 0.f);

        if (m_profileManager.GetRecoStartDirectionAndPosition(evt, pfparticle, recoStartDirection, recoPosition))
        {
            ProfileManager::EnergyProfiles energyProfiles_recoStart = m_profileManager.CreateProfiles(evt, pfparticle, recoStartDirection, recoPosition);
            m_longitudinal_recoStart = energyProfiles_recoStart.m_longitudinal;
            m_transverse_recoStart = energyProfiles_recoStart.m_transverse;
            m_energy_recoStart = energyProfiles_recoStart.m_energy;
            m_fitDirection_recoStart = std::vector<float>({static_cast<float>(recoStartDirection.X()), static_cast<float>(recoStartDirection.Y()), static_cast<float>(recoStartDirection.Z())});

            m_profileManager.FillFitPositions(evt, pfparticle, recoStartDirection, recoPosition, m_fitPositions_recoStart);
        }

        // PCA method
        TVector3 pcaPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 pcaDirection = TVector3(0.f, 0.f, 0.f);

        if (m_profileManager.GetPCADirectionAndPosition(evt, pfparticle, pcaDirection, pcaPosition))
        {
            ProfileManager::EnergyProfiles energyProfiles_pca = m_profileManager.CreateProfiles(evt, pfparticle, pcaDirection, pcaPosition);
            m_longitudinal_pca = energyProfiles_pca.m_longitudinal;
            m_transverse_pca = energyProfiles_pca.m_transverse;
            m_energy_pca = energyProfiles_pca.m_energy;
            m_fitDirection_pca = std::vector<float>({static_cast<float>(pcaDirection.X()), static_cast<float>(pcaDirection.Y()), static_cast<float>(pcaDirection.Z())});

            m_profileManager.FillFitPositions(evt, pfparticle, pcaDirection, pcaPosition, m_fitPositions_pca);
        }

        // Initial track stub method
        TVector3 trackStubPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 trackStubDirection = TVector3(0.f, 0.f, 0.f);

        if (m_profileManager.GetTrackStubDirectionAndPosition(evt, pfparticle, trackStubDirection, trackStubPosition))
        {
            ProfileManager::EnergyProfiles energyProfiles_trackStub = m_profileManager.CreateProfiles(evt, pfparticle, trackStubDirection, trackStubPosition);
            m_longitudinal_trackStub = energyProfiles_trackStub.m_longitudinal;
            m_transverse_trackStub = energyProfiles_trackStub.m_transverse;
            m_energy_trackStub = energyProfiles_trackStub.m_energy;
            m_fitDirection_trackStub = std::vector<float>({static_cast<float>(trackStubDirection.X()), static_cast<float>(trackStubDirection.Y()), static_cast<float>(trackStubDirection.Z())});

            m_profileManager.FillFitPositions(evt, pfparticle, trackStubDirection, trackStubPosition, m_fitPositions_trackStub);
        }

        // Sliding linear fit (essentially a track fit?)
        ProfileManager::EnergyProfiles energyProfiles_track = m_profileManager.CreateProfilesFromTrack(evt, pfparticle);
        m_longitudinal_track = energyProfiles_track.m_longitudinal;
        m_transverse_track = energyProfiles_track.m_transverse;
        m_energy_track = energyProfiles_track.m_energy;

        m_profileManager.FillFitPositionsFromTrack(evt, pfparticle, m_fitPositions_track);

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

void VenusaurusTrainingFiles::Reset()
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

  m_longitudinal_true.clear();
  m_transverse_true.clear();
  m_energy_true.clear();
  m_fitDirection_true.clear();
  m_fitPositions_true.clear();

  m_longitudinal_recoStart.clear();
  m_transverse_recoStart.clear();
  m_energy_recoStart.clear(); 
  m_fitDirection_recoStart.clear();
  m_fitPositions_recoStart.clear();

  m_longitudinal_pca.clear();
  m_transverse_pca.clear();
  m_energy_pca.clear(); 
  m_fitDirection_pca.clear();
  m_fitPositions_pca.clear();

  m_longitudinal_trackStub.clear();
  m_transverse_trackStub.clear();
  m_energy_trackStub.clear(); 
  m_fitDirection_trackStub.clear();
  m_fitPositions_trackStub.clear();

  m_longitudinal_track.clear();
  m_transverse_track.clear();
  m_energy_track.clear(); 
  m_fitPositions_track.clear();

  m_longitudinalProfile.clear();
  m_transverseProfile.clear();

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

void VenusaurusTrainingFiles::beginJob()
{
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    m_tree = tfs->make<TTree>("venusaur", "venusaur");

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

    m_tree->Branch("Longitudinal_true", &m_longitudinal_true);
    m_tree->Branch("Transverse_true", &m_transverse_true);
    m_tree->Branch("Energy_true", &m_energy_true);
    m_tree->Branch("FitDirection_true", &m_fitDirection_true);
    m_tree->Branch("FitPositions_true", &m_fitPositions_true);

    m_tree->Branch("Longitudinal_recoStart", &m_longitudinal_recoStart);
    m_tree->Branch("Transverse_recoStart", &m_transverse_recoStart);
    m_tree->Branch("Energy_recoStart", &m_energy_recoStart);
    m_tree->Branch("FitDirection_recoStart", &m_fitDirection_recoStart);
    m_tree->Branch("FitPositions_recoStart", &m_fitPositions_recoStart);

    m_tree->Branch("Longitudinal_pca", &m_longitudinal_pca);
    m_tree->Branch("Transverse_pca", &m_transverse_pca);
    m_tree->Branch("Energy_pca", &m_energy_pca);
    m_tree->Branch("FitDirection_pca", &m_fitDirection_pca);
    m_tree->Branch("FitPositions_pca", &m_fitPositions_pca);

    m_tree->Branch("Longitudinal_trackStub", &m_longitudinal_trackStub);
    m_tree->Branch("Transverse_trackStub", &m_transverse_trackStub);
    m_tree->Branch("Energy_trackStub", &m_energy_trackStub);
    m_tree->Branch("FitDirection_trackStub", &m_fitDirection_trackStub);
    m_tree->Branch("FitPositions_trackStub", &m_fitPositions_trackStub);

    m_tree->Branch("Longitudinal_track", &m_longitudinal_track);
    m_tree->Branch("Transverse_track", &m_transverse_track);
    m_tree->Branch("Energy_track", &m_energy_track);
    m_tree->Branch("FitPositions_track", &m_fitPositions_track);

    m_tree->Branch("LongitudinalProfile", &m_longitudinalProfile);
    m_tree->Branch("TransverseProfile", &m_transverseProfile);

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

void VenusaurusTrainingFiles::endJob()
{
}


} //namespace ivysaurus

