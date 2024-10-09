///////////////////////////////////////////////
// CCNuSelection analyzer
//
// Creates tree on which to perform Pandora-based numu/nue selection
// I Mawby, D Brailsford May 2023
///////////////////////////////////////////////

#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//STL
#include <iostream>

//ROOT
#include "TTree.h"
#include "TMVA/Reader.h"

//ART
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

//LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

//DUNE
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/TrackPID/algorithms/CTPHelper.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

//Custom
#include "dunereco/FDSelections/EnergyReco/TrueEnergyCalc.h"
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -999;
constexpr double kDefDoub = -999.0;
constexpr int kMaxMCParticles = 1000;
constexpr int kMaxPFParticles = 100;
constexpr int kMaxDEDXPoints = 1000;

namespace FDSelection {
  class CCNuSelection;
}

class FDSelection::CCNuSelection : public art::EDAnalyzer {
public:
  explicit CCNuSelection(fhicl::ParameterSet const & p);
  CCNuSelection(CCNuSelection const &) = delete;
  CCNuSelection(CCNuSelection &&) = delete;
  CCNuSelection & operator = (CCNuSelection const &) = delete;
  CCNuSelection & operator = (CCNuSelection &&) = delete;
  void analyze(art::Event const & e) override;
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr) override;
  void endSubRun(art::SubRun const & sr) override;
  void endJob() override;

private:
  void Reset();
  void FillPandoraMaps(art::Event const& evt);
  void GetTruthInfo(art::Event const & evt);
  void FillMCParticleInfo(art::Event const & evt);
  void FillVertexInfo(art::Event const & evt);
  void FillPFParticleInfo(art::Event const & evt);
  void FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);
  void FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);

  ////////////////////////////////////////
  // Pandora maps
  ////////////////////////////////////////
  lar_pandora::MCParticleMap fMCParticleMap;   // Linking TrackID -> MCParticle
  lar_pandora::PFParticleMap fPFPMap;          // Linking Self() -> PFParticle
  ////////////////////////////////////////
  // Trees
  ////////////////////////////////////////
  TTree *fTree;
  ////////////////////////////////////////
  // Event info
  ////////////////////////////////////////
  // Event identification
  int fRun;
  int fSubRun;
  int fEvent;
  // Neutrino 
  int fNuPdg;        // Interaction PDG
  int fBeamPdg;      // PDG at point of creation
  int fNuTrackID;
  int fNC;           // 1=is NC, 0=otherwise
  int fMode;         // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fTargetZ;      // Atomic number of scattering target
  double fENu; 
  double fNuMomX;    // Neutrino momentum
  double fNuMomY;
  double fNuMomZ;
  double fNuX;       // Interaction positions
  double fNuY;
  double fNuZ;
  ////////////////////////////////////////
  // Event-level reco
  ////////////////////////////////////////
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
  int fNMCParticles;
  int fNRecoPFPs;
  ////////////////////////////////////////
  // MCParticle Info  
  ////////////////////////////////////////
  int fMCTrackID[kMaxMCParticles];
  int fMCParentTrackID[kMaxMCParticles];
  ////////////////////////////////////////
  // PFParticle Info  
  ////////////////////////////////////////
  // Truth
  int fRecoPFPTruePDG[kMaxPFParticles];
  int fRecoPFPTrueTrackID[kMaxPFParticles];
  bool fRecoPFPTruePrimary[kMaxPFParticles];
  int fRecoPFPTrueGeneration[kMaxPFParticles];
  int fRecoPFPTrueParentTrackID[kMaxPFParticles];
  int fRecoPFPTrueParentPDG[kMaxPFParticles];
  int fRecoPFPTrueVisibleGeneration[kMaxPFParticles];
  int fRecoPFPTrueVisibleParentTrackID[kMaxPFParticles];
  int fRecoPFPTrueVisibleParentPDG[kMaxPFParticles];
  int fRecoPFPTrueVisibleParentSelf[kMaxPFParticles];
  int fRecoPFPTrueVisibleParentPFPIndex[kMaxPFParticles];
  // Reco
  int fRecoPFPSelf[kMaxPFParticles];
  bool fRecoPFPIsPrimary[kMaxPFParticles];
  int fRecoPFPRecoGeneration[kMaxPFParticles];
  int fRecoPFPRecoParentSelf[kMaxPFParticles];
  double fRecoPFPTrackShowerScore[kMaxPFParticles];
  double fRecoPFPVertexX[kMaxPFParticles];
  double fRecoPFPVertexY[kMaxPFParticles];
  double fRecoPFPVertexZ[kMaxPFParticles];
  int fRecoPFPRecoNHits[kMaxPFParticles];
  int fRecoPFPRecoNSpacepoints[kMaxPFParticles];
  std::vector<std::vector<double>> fRecoPFPSpacepointX;
  std::vector<std::vector<double>> fRecoPFPSpacepointY;
  std::vector<std::vector<double>> fRecoPFPSpacepointZ;
  double fRecoPFPRecoCompleteness[kMaxPFParticles];
  double fRecoPFPRecoHitPurity[kMaxPFParticles];
  bool fRecoPFPTrackFitSuccess[kMaxPFParticles];
  bool fRecoPFPShowerFitSuccess[kMaxPFParticles];
  ////////////////////////////////////////
  // Track info
  ////////////////////////////////////////
  // Reco
  double fRecoTrackRecoStartX[kMaxPFParticles];
  double fRecoTrackRecoStartY[kMaxPFParticles];
  double fRecoTrackRecoStartZ[kMaxPFParticles];
  double fRecoTrackRecoEndX[kMaxPFParticles];
  double fRecoTrackRecoEndY[kMaxPFParticles];
  double fRecoTrackRecoEndZ[kMaxPFParticles];
  double fRecoTrackRecoLength[kMaxPFParticles];
  // Pandizzle PID (minus michel vars...)
  double fRecoTrackDeflecAngleSD[kMaxPFParticles];
  double fRecoTrackLength[kMaxPFParticles];
  double fRecoTrackEvalRatio[kMaxPFParticles];
  double fRecoTrackConcentration[kMaxPFParticles];
  double fRecoTrackCoreHaloRatio[kMaxPFParticles];
  double fRecoTrackConicalness[kMaxPFParticles];
  double fRecoTrackdEdxStart[kMaxPFParticles];
  double fRecoTrackdEdxEnd[kMaxPFParticles];
  double fRecoTrackdEdxEndRatio[kMaxPFParticles];
  int fNTrajPoints;
  int fTrajPointStartIndex[kMaxPFParticles];
  int fTrajPointEndIndex[kMaxPFParticles];
  double fRecoTrackRecodEdx[kMaxDEDXPoints];
  double fRecoTrackRecoRR[kMaxDEDXPoints];
  ////////////////////////////////////////
  // Shower Info
  ////////////////////////////////////////
  double fRecoShowerRecoStartX[kMaxPFParticles];
  double fRecoShowerRecoStartY[kMaxPFParticles];
  double fRecoShowerRecoStartZ[kMaxPFParticles];
  double fRecoShowerRecoLength[kMaxPFParticles];
  double fRecoShowerRecoOpeningAngle[kMaxPFParticles];
  double fRecoShowerRecodEdx[kMaxPFParticles][3];
  int fRecoShowerRecoBestPlane[kMaxPFParticles];
  double fRecoShowerRecoEnergy[kMaxPFParticles][3];
  ////////////////////////////////////////
  //Module labels
  ////////////////////////////////////////
  std::string fNuGenModuleLabel;
  std::string fLargeantModuleLabel;
  std::string fWireModuleLabel;
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fRecoModuleLabel;
  std::string fHitsModuleLabel;
  ////////////////////////////////////////
  //Algs
  ////////////////////////////////////////
  PandizzleAlg fPandizzleAlg;
  PandrizzleAlg fPandrizzleAlg;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
};

//////////////////////////////////////////////////////////////////////////////////

FDSelection::CCNuSelection::CCNuSelection(fhicl::ParameterSet const & pset) :
  EDAnalyzer(pset),
  fNuGenModuleLabel(pset.get<std::string>("NuGenModuleLabel")),
  fLargeantModuleLabel(pset.get<std::string>("LargeantModuleLabel")),
  fWireModuleLabel(pset.get<std::string>("WireModuleLabel")),
  fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
  fRecoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
  fHitsModuleLabel(pset.get<std::string>("HitsModuleLabel")),
  fPandizzleAlg(pset.get<fhicl::ParameterSet>("PandizzleConfig")),
  fPandrizzleAlg(pset.get<fhicl::ParameterSet>("PandrizzleConfig")),
  fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"), fTrackModuleLabel, fShowerModuleLabel,
      fHitsModuleLabel, fWireModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fRecoModuleLabel)
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
    //std::cout << "AAAAAA" << std::endl;
    Reset();
    //std::cout << "BBB" << std::endl;
    fRun = evt.run();
    fSubRun = evt.subRun();
    fEvent = evt.event();
    //std::cout << "CCC" << std::endl;
    FillPandoraMaps(evt);
    //std::cout << "EEE" << std::endl;
    GetTruthInfo(evt);
    //std::cout << "FFF" << std::endl;
    FillVertexInfo(evt);
    //std::cout << "GGG" << std::endl;
    FillPFParticleInfo(evt);
    //std::cout << "HHH" << std::endl;

    fTree->Fill();
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::beginJob()
{
    ///////////////////////////
    // POT tree
    ///////////////////////////
    art::ServiceHandle<art::TFileService> tfs;

    ///////////////////////////
    // CCNuSelection tree
    ///////////////////////////
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");

    ////////////////////////////
    // Event info
    ////////////////////////////
    fTree->Branch("Event_Run", &fRun);
    fTree->Branch("Event_Subrun", &fSubRun);
    fTree->Branch("Event_Event", &fEvent);
    ////////////////////////////
    // True info
    ////////////////////////////
    fTree->Branch("Nu_True_PDG", &fNuPdg);
    fTree->Branch("Nu_True_SimID", &fNuTrackID);
    fTree->Branch("Nu_True_IsNC", &fNC);
    fTree->Branch("Nu_True_Mode", &fMode);
    fTree->Branch("Nu_True_TargetZ", &fTargetZ);
    fTree->Branch("Nu_True_Energy", &fENu);
    fTree->Branch("Nu_True_MomX", &fNuMomX);
    fTree->Branch("Nu_True_MomY", &fNuMomY);
    fTree->Branch("Nu_True_MomZ", &fNuMomZ);
    fTree->Branch("Nu_True_VertexX", &fNuX);
    fTree->Branch("Nu_True_VertexY", &fNuY);
    fTree->Branch("Nu_True_VertexZ", &fNuZ);
    ////////////////////////////
    // Event-level reco info 
    ////////////////////////////
    fTree->Branch("Event_NMCParticles", &fNMCParticles);
    fTree->Branch("Event_NRecoPFPs", &fNRecoPFPs);
    fTree->Branch("Nu_Reco_VertexX", &fRecoNuVtxX);
    fTree->Branch("Nu_Reco_VertexY", &fRecoNuVtxY);
    fTree->Branch("Nu_Reco_VertexZ", &fRecoNuVtxZ);
    ////////////////////////////
    // MCParticle info
    ////////////////////////////
    fTree->Branch("MC_SimID", fMCTrackID, "MC_SimID[Event_NMCParticles]/I");
    fTree->Branch("MC_ParentSimID", fMCParentTrackID, "MC_ParentSimID[Event_NMCParticles]/I");
    ////////////////////////////
    // PFParticle info
    ////////////////////////////
    // ID
    fTree->Branch("PFP_RecoID", fRecoPFPSelf, "PFP_RecoID[Event_NRecoPFPs]/I");
    // True
    fTree->Branch("MatchedMC_PDG", fRecoPFPTruePDG, "MatchedMC_PDG[Event_NRecoPFPs]/I");
    fTree->Branch("MatchedMC_SimID", fRecoPFPTrueTrackID, "MatchedMC_SimID[Event_NRecoPFPs]/I");
    fTree->Branch("MatchedMC_IsTruePrimary", fRecoPFPTruePrimary,"MatchedMC_IsTruePrimary[Event_NRecoPFPs]/O");
    fTree->Branch("MatchedMC_Generation", fRecoPFPTrueGeneration, "MatchedMC_Generation[Event_NRecoPFPs]/I");
    fTree->Branch("MatchedMC_ParentSimID", fRecoPFPTrueParentTrackID, "MatchedMC_ParentSimID[Event_NRecoPFPs]/I");
    fTree->Branch("MatchedMC_ParentPDG", fRecoPFPTrueParentPDG, "MatchedMC_ParentPDG[Event_NRecoPFPs]/I");
    // Reco
    // Current Pandora hierarchy stuff
    fTree->Branch("PFP_CurrentReco_IsPrimary", fRecoPFPIsPrimary, "PFP_CurrentReco_IsPrimary[Event_NRecoPFPs]/O");
    fTree->Branch("PFP_CurrentReco_Generation", fRecoPFPRecoGeneration, "PFP_CurrentReco_Generation[Event_NRecoPFPs]/I");
    fTree->Branch("PFP_CurrentReco_ParentRecoID", fRecoPFPRecoParentSelf, "PFP_CurrentReco_ParentRecoID[Event_NRecoPFPs]/I");
    fTree->Branch("PFP_TrackShowerScore", fRecoPFPTrackShowerScore, "PFP_TrackShowerScore[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_VertexX", fRecoPFPVertexX, "PFP_VertexX[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_VertexY", fRecoPFPVertexY, "PFP_VertexY[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_VertexZ", fRecoPFPVertexZ, "PFP_VertexZ[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_N2DHits", fRecoPFPRecoNHits,"PFP_N2DHits[Event_NRecoPFPs]/I");
    fTree->Branch("PFP_NSpacepoints", fRecoPFPRecoNSpacepoints, "PFP_NSpacepoints[Event_NRecoPFPs]/I");
    fTree->Branch("PFP_SpacepointX", &fRecoPFPSpacepointX);
    fTree->Branch("PFP_SpacepointY", &fRecoPFPSpacepointY);
    fTree->Branch("PFP_SpacepointZ", &fRecoPFPSpacepointZ);
    fTree->Branch("PFP_Completeness", fRecoPFPRecoCompleteness, "PFP_Completeness[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_Purity", fRecoPFPRecoHitPurity, "PFP_Purity[Event_NRecoPFPs]/D");
    fTree->Branch("PFP_ShowerFitSuccess", fRecoPFPShowerFitSuccess, "PFP_ShowerFitSuccess[Event_NRecoPFPs]/O");
    fTree->Branch("PFP_TrackFitSuccess", fRecoPFPTrackFitSuccess, "PFP_TrackFitSuccess[Event_NRecoPFPs]/O");
    ////////////////////////////
    // Track Info
    ////////////////////////////
    // Reco
    fTree->Branch("Track_StartX", fRecoTrackRecoStartX, "Track_StartX[Event_NRecoPFPs]/D");
    fTree->Branch("Track_StartY", fRecoTrackRecoStartY, "Track_StartY[Event_NRecoPFPs]/D");
    fTree->Branch("Track_StartZ", fRecoTrackRecoStartZ, "Track_StartZ[Event_NRecoPFPs]/D");
    fTree->Branch("Track_EndX", fRecoTrackRecoEndX, "Track_EndX[Event_NRecoPFPs]/D");
    fTree->Branch("Track_EndY", fRecoTrackRecoEndY, "Track_EndY[Event_NRecoPFPs]/D");
    fTree->Branch("Track_EndZ", fRecoTrackRecoEndZ, "Track_EndZ[Event_NRecoPFPs]/D");
    fTree->Branch("Track_Length", fRecoTrackRecoLength, "Track_Length[Event_NRecoPFPs]/D");
    // Pandizzle
    fTree->Branch("Track_DeflecAngleSD", fRecoTrackDeflecAngleSD, "Track_DeflecAngleSD[Event_NRecoPFPs]/D");
    fTree->Branch("Track_Length", fRecoTrackLength, "Track_Length[Event_NRecoPFPs]/D");
    fTree->Branch("Track_EvalRatio", fRecoTrackEvalRatio, "Track_EvalRatio[Event_NRecoPFPs]/D");
    fTree->Branch("Track_Concentration", fRecoTrackConcentration, "Track_Concentration[Event_NRecoPFPs]/D");
    fTree->Branch("Track_CoreHaloRatio", fRecoTrackCoreHaloRatio, "Track_CoreHaloRatio[Event_NRecoPFPs]/D");
    fTree->Branch("Track_Conicalness", fRecoTrackConicalness, "Track_Conicalness[Event_NRecoPFPs]/D");
    fTree->Branch("Track_dEdxStart", fRecoTrackdEdxStart, "Track_dEdxStart[Event_NRecoPFPs]/D");
    fTree->Branch("Track_dEdxEnd", fRecoTrackdEdxEnd, "Track_dEdxEnd[Event_NRecoPFPs]/D");
    fTree->Branch("Track_dEdxEndRatio", fRecoTrackdEdxEndRatio, "Track_dEdxEndRatio[Event_NRecoPFPs]/D");
    fTree->Branch("Track_NTrajPoints", &fNTrajPoints);
    fTree->Branch("Track_TrajPointStartIndex", fTrajPointStartIndex, "Track_TrajPointStartIndex[Event_NRecoPFPs]/I");
    fTree->Branch("Track_TrajPointEndIndex", fTrajPointEndIndex, "Track_TrajPointEndIndex[Event_NRecoPFPs]/I");
    fTree->Branch("Track_dEdX", fRecoTrackRecodEdx, "Track_dEdX[Track_NTrajPoints]/D");
    fTree->Branch("Track_RR", fRecoTrackRecoRR, "Track_RR[Track_NTrajPoints]/D");
    ///////////////////////////
    // Shower Info
    ///////////////////////////
    // Reco
    fTree->Branch("Shower_StartX", fRecoShowerRecoStartX, "Shower_StartX[Event_NRecoPFPs]/D");
    fTree->Branch("Shower_StartY", fRecoShowerRecoStartY, "Shower_StartY[Event_NRecoPFPs]/D");
    fTree->Branch("Shower_StartZ", fRecoShowerRecoStartZ, "Shower_StartZ[Event_NRecoPFPs]/D");
    fTree->Branch("Shower_Length", &fRecoShowerRecoLength, "Shower_Length[Event_NRecoPFPs]/D");
    fTree->Branch("Shower_OpeningAngle", &fRecoShowerRecoOpeningAngle, "Shower_OpeningAngle[Event_NRecoPFPs]/D");
    fTree->Branch("Shower_InitialdEdx", fRecoShowerRecodEdx, "Shower_InitialdEdx[Event_NRecoPFPs][3]/D");
    fTree->Branch("Shower_BestPlane", &fRecoShowerRecoBestPlane, "Shower_BestPlane[Event_NRecoPFPs]/I");
    fTree->Branch("Shower_Energy", fRecoShowerRecoEnergy, "Shower_Energy[Event_NRecoPFPs][3]/D");
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::beginSubRun(art::SubRun const & sr)
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::endSubRun(const art::SubRun& sr)
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::endJob()
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::Reset()
{
    ////////////////////////////////////////
    // Pandora maps
    ////////////////////////////////////////
    fMCParticleMap.clear();
    fPFPMap.clear();
    ////////////////////////////
    // Event info
    ////////////////////////////
    fRun = kDefInt;
    fSubRun = kDefInt;
    fEvent = kDefInt;
    ////////////////////////////
    // True info
    ////////////////////////////
    fNuPdg = kDefInt; 
    fNuTrackID = kDefInt;
    fNC = kDefInt;    
    fMode = kDefInt; 
    fTargetZ = kDefInt;
    fENu = kDefDoub; 
    fNuMomX = kDefDoub; 
    fNuMomY = kDefDoub;
    fNuMomZ = kDefDoub;
    fNuX = kDefDoub; 
    fNuY = kDefDoub;
    fNuZ = kDefDoub;
    ////////////////////////////
    // Event-level reco info 
    ////////////////////////////
    fNMCParticles = 0;
    fNRecoPFPs = 0;
    fRecoNuVtxX = kDefDoub;
    fRecoNuVtxY = kDefDoub;
    fRecoNuVtxZ = kDefDoub;

    fRecoPFPSpacepointX.clear();
    fRecoPFPSpacepointY.clear();
    fRecoPFPSpacepointZ.clear();

    fRecoPFPSpacepointX.resize(0);
    fRecoPFPSpacepointY.resize(0);
    fRecoPFPSpacepointZ.resize(0);

    ////////////////////////////
    // MCParticle info
    ////////////////////////////
    for (int i = 0; i < kMaxMCParticles; i++)
    {
        fMCTrackID[i] = kDefDoub;
        fMCParentTrackID[i] = kDefDoub;
    }
    ////////////////////////////
    // PFParticle info
    ////////////////////////////
    for (int i = 0; i < kMaxPFParticles; i++)
    {
        // ID
        fRecoPFPSelf[i] = kDefInt;
        // Truth
        fRecoPFPTruePDG[i] = kDefInt;
        fRecoPFPTrueTrackID[i] = kDefInt;
        fRecoPFPTruePrimary[i] = false;
        fRecoPFPTrueGeneration[i] = kDefInt;
        fRecoPFPTrueParentTrackID[i] = kDefInt;
        fRecoPFPTrueParentPDG[i] = kDefInt;
        // Reco
        fRecoPFPIsPrimary[i] = false;
        fRecoPFPRecoGeneration[i] = kDefInt;
        fRecoPFPRecoParentSelf[i] = kDefInt;
        fRecoPFPTrackShowerScore[i] = kDefDoub;
        fRecoPFPVertexX[i] = kDefDoub;
        fRecoPFPVertexY[i] = kDefDoub;
        fRecoPFPVertexZ[i] = kDefDoub;
        fRecoPFPRecoNHits[i] = kDefInt;
        fRecoPFPRecoNSpacepoints[i] = kDefInt;
        fRecoPFPShowerFitSuccess[i] = false;
        fRecoPFPTrackFitSuccess[i] = false;
        fRecoPFPRecoCompleteness[i] = kDefDoub;
        fRecoPFPRecoHitPurity[i] = kDefDoub;
        ////////////////////////////
        // Track stuff
        ////////////////////////////
        // Reco
        fRecoTrackRecoStartX[i] = kDefDoub;
        fRecoTrackRecoStartY[i] = kDefDoub;
        fRecoTrackRecoStartZ[i] = kDefDoub;
        fRecoTrackRecoEndX[i] = kDefDoub;
        fRecoTrackRecoEndY[i] = kDefDoub;
        fRecoTrackRecoEndZ[i] = kDefDoub;
        fRecoTrackRecoLength[i] = kDefDoub;
        // Pandizzle
        fRecoTrackDeflecAngleSD[i] = kDefDoub;
        fRecoTrackLength[i] = kDefDoub;
        fRecoTrackEvalRatio[i] = kDefDoub;
        fRecoTrackConcentration[i] = kDefDoub;
        fRecoTrackCoreHaloRatio[i] = kDefDoub;
        fRecoTrackConicalness[i] = kDefDoub;
        fRecoTrackdEdxStart[i] = kDefDoub;
        fRecoTrackdEdxEnd[i] = kDefDoub;
        fRecoTrackdEdxEndRatio[i] = kDefDoub;

        fNTrajPoints = 0;
        fTrajPointStartIndex[i] = kDefInt;
        fTrajPointEndIndex[i] = kDefInt;

        for (int j = 0; j < kMaxDEDXPoints; ++j)
        {
            fRecoTrackRecodEdx[j] = kDefDoub;
            fRecoTrackRecoRR[j] = kDefDoub;
        }

        ////////////////////////////
        // Shower stuff
        ////////////////////////////
        // Reco
        fRecoShowerRecoStartX[i] = kDefDoub;
        fRecoShowerRecoStartY[i] = kDefDoub;
        fRecoShowerRecoStartZ[i] = kDefDoub;
        fRecoShowerRecoLength[i] = kDefDoub;
        fRecoShowerRecoOpeningAngle[i] = kDefDoub;
        fRecoShowerRecoBestPlane[i] = kDefInt;
        for (int j = 0; j < 3; j++)
        {
            fRecoShowerRecodEdx[i][j] = kDefDoub;
            fRecoShowerRecoEnergy[i][j] = kDefDoub;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    const std::vector<art::Ptr<simb::MCParticle>> mcParticles = dune_ana::DUNEAnaEventUtils::GetMCParticles(evt, fLargeantModuleLabel);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticles, fMCParticleMap);

    // PFParticle map
    std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfps, fPFPMap);
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetTruthInfo(art::Event const & evt)
{
    const std::vector<art::Ptr<simb::MCTruth>> mcTruths = dune_ana::DUNEAnaEventUtils::GetMCTruths(evt, fNuGenModuleLabel);

    if (mcTruths.empty())
        return;

    art::Ptr<simb::MCTruth> mcTruth = mcTruths.at(0);

    if (mcTruth->Origin() != simb::kBeamNeutrino)
        return;

    // Neutrino
    const simb::MCNeutrino &mcNeutrino = mcTruth->GetNeutrino();
    fNuTrackID = mcNeutrino.Nu().TrackId();
    fNuPdg = mcNeutrino.Nu().PdgCode();
    fNC = mcNeutrino.CCNC();
    fMode = mcNeutrino.Mode();
    fTargetZ = mcNeutrino.Target()%100000000/10000;
    fENu = mcNeutrino.Nu().E();
    fNuX = mcNeutrino.Nu().Vx();
    fNuY = mcNeutrino.Nu().Vy();
    fNuZ = mcNeutrino.Nu().Vz();
    fNuMomX = mcNeutrino.Nu().Momentum().X();
    fNuMomY = mcNeutrino.Nu().Momentum().Y();
    fNuMomZ = mcNeutrino.Nu().Momentum().Z();

    // Fill MCParticle vectors
    for (auto &entry : fMCParticleMap)
    {
        if (fNMCParticles == kMaxMCParticles)
            break;

        fMCTrackID[fNMCParticles] = entry.first;
        fMCParentTrackID[fNMCParticles] = entry.second->Mother();

        ++fNMCParticles;
    }
}

///////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillVertexInfo(art::Event const & evt)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
        return;

    art::Ptr<recob::PFParticle> nu_pfp = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fRecoModuleLabel);

    try
    {
        art::Ptr<recob::Vertex> nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nu_pfp, evt, fRecoModuleLabel);

        fRecoNuVtxX = nuVertex->position().X();
        fRecoNuVtxY = nuVertex->position().Y();
        fRecoNuVtxZ = nuVertex->position().Z();
    }
    catch(...)
    {
    }

    return;
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillPFParticleInfo(art::Event const & evt)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
        return;

    std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);
    std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);

    int pfpIndex = -1;
    fNRecoPFPs = 0;

    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        // Skip the neutrino
        if ((std::fabs(pfp->PdgCode()) == 12) || (std::fabs(pfp->PdgCode()) == 14) || (std::fabs(pfp->PdgCode()) == 16))
            continue;

        pfpIndex++;
        fNRecoPFPs++;

        if (pfpIndex == kMaxPFParticles)
            break;

        // Simple things
        fRecoPFPSelf[pfpIndex] = pfp->Self();
        fRecoPFPRecoGeneration[pfpIndex] = lar_pandora::LArPandoraHelper::GetGeneration(fPFPMap, pfp);
        fRecoPFPIsPrimary[pfpIndex] = (fRecoPFPRecoGeneration[pfpIndex] == 2 ? true : false);

        // Vertex
        try
        {
            art::Ptr<recob::Vertex> recoVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, fRecoModuleLabel);

            fRecoPFPVertexX[pfpIndex] = recoVertex->position().X();
            fRecoPFPVertexY[pfpIndex] = recoVertex->position().Y();
            fRecoPFPVertexZ[pfpIndex] = recoVertex->position().Z();
        }
        catch(...)
        {
        }

        // Track/shower score
        try
        {
            const art::Ptr<larpandoraobj::PFParticleMetadata> metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfp, evt, fRecoModuleLabel);

            if (metadata->GetPropertiesMap().find("TrackScore") != metadata->GetPropertiesMap().end())
                fRecoPFPTrackShowerScore[pfpIndex] = metadata->GetPropertiesMap().at("TrackScore");
        }
        catch (...) {}

        // 2D hits
        const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fRecoModuleLabel);
        fRecoPFPRecoNHits[pfpIndex] = pfpHits.size();

        // 3D spacepoints
        const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, fRecoModuleLabel);
        fRecoPFPRecoNSpacepoints[pfpIndex] = spacepoints.size();

        std::vector<double> spacepointX, spacepointY, spacepointZ;

        for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
        {
            spacepointX.push_back(spacepoint->XYZ()[0]);
            spacepointY.push_back(spacepoint->XYZ()[1]);
            spacepointZ.push_back(spacepoint->XYZ()[2]);
        }

        fRecoPFPSpacepointX.push_back(spacepointX);
        fRecoPFPSpacepointY.push_back(spacepointY);
        fRecoPFPSpacepointZ.push_back(spacepointZ);

        // Truth information
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
        int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);

        if (TruthMatchUtils::Valid(g4id))
        {
            fRecoPFPRecoCompleteness[pfpIndex] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, pfpHits, eventHitList, g4id);
            fRecoPFPRecoHitPurity[pfpIndex] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, pfpHits, g4id);

            art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
            const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);

            if (matched_mcparticle)
            {
                fRecoPFPTruePDG[pfpIndex] = matched_mcparticle->PdgCode();
                fRecoPFPTrueTrackID[pfpIndex] = matched_mcparticle->TrackId();

                if (matched_mcparticle->Mother() == 0) 
                    fRecoPFPTruePrimary[pfpIndex] = true;
                else 
                    fRecoPFPTruePrimary[pfpIndex] = false;
            }
        }

        // Fill the track & shower information
        FillRecoTrackInfo(evt, pfp, pfpIndex);
        FillRecoShowerInfo(evt, pfp, pfpIndex);
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpIndex)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fRecoModuleLabel, fTrackModuleLabel))
        return;

    fRecoPFPTrackFitSuccess[pfpIndex] = true;

    art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, fRecoModuleLabel, fTrackModuleLabel);

    // General track variables
    recob::Track::Point_t trackStart, trackEnd;
    std::tie(trackStart, trackEnd) = track->Extent(); 
    fRecoTrackRecoStartX[pfpIndex] = trackStart.X();
    fRecoTrackRecoStartY[pfpIndex] = trackStart.Y();
    fRecoTrackRecoStartZ[pfpIndex] = trackStart.Z();
    fRecoTrackRecoEndX[pfpIndex] = trackEnd.X();
    fRecoTrackRecoEndY[pfpIndex] = trackEnd.Y();
    fRecoTrackRecoEndZ[pfpIndex] = trackEnd.Z();
    fRecoTrackRecoLength[pfpIndex] = track->Length();

    // Pandizzle variables
    FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(track, evt));
    fRecoTrackDeflecAngleSD[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
    fRecoTrackLength[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
    fRecoTrackEvalRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
    fRecoTrackConcentration[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
    fRecoTrackCoreHaloRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
    fRecoTrackConicalness[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
    fRecoTrackdEdxStart[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
    fRecoTrackdEdxEnd[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
    fRecoTrackdEdxEndRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);

    // Primary Track Calorimetry
    try
    {
        const art::Ptr<anab::Calorimetry> &trackCalo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(track, evt, fTrackModuleLabel, "pandoracalo");

        int nTrajPoints = trackCalo->dEdx().size();
        int trajPointStartIndex = fNTrajPoints;
        int trajPointEndIndex = trajPointStartIndex + (nTrajPoints - 1);

        if (trajPointEndIndex < kMaxDEDXPoints)
        {
            fNTrajPoints += nTrajPoints;
            fTrajPointStartIndex[pfpIndex] = trajPointStartIndex;
            fTrajPointEndIndex[pfpIndex] = trajPointEndIndex;

            for (int i = 0; i < nTrajPoints; i++)
            {
                fRecoTrackRecodEdx[trajPointStartIndex + i] = trackCalo->dEdx()[i];
                fRecoTrackRecoRR[trajPointStartIndex + i] = trackCalo->ResidualRange()[i];
            }
        }
    }
    catch(...)
    {
    }
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpIndex)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel))
        return;

    fRecoPFPShowerFitSuccess[pfpIndex] = true;

    art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel);

    // General shower variables
    fRecoShowerRecoStartX[pfpIndex] = shower->ShowerStart().X();
    fRecoShowerRecoStartY[pfpIndex] = shower->ShowerStart().Y();
    fRecoShowerRecoStartZ[pfpIndex] = shower->ShowerStart().Z();
    fRecoShowerRecoBestPlane[pfpIndex] = shower->best_plane();
    fRecoShowerRecoLength[pfpIndex] = shower->Length();
    fRecoShowerRecoOpeningAngle[pfpIndex] = shower->OpenAngle();

    // Momentum and energy
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(shower, evt)));

    if (shower->dEdx().size() > 0)
    {
      for (int i_plane = 0; i_plane < 3; i_plane++)
      {
        fRecoShowerRecodEdx[pfpIndex][i_plane] = shower->dEdx()[i_plane];
        fRecoShowerRecoEnergy[pfpIndex][i_plane] = shower->Energy()[i_plane];
      }
    }
}

//////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(FDSelection::CCNuSelection)
