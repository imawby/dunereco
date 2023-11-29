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
#include "dunereco/TrackPID/products/CTPResult.h"

#include "dunereco/CVN/func/Result.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

//Custom
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);
constexpr int kDefMaxNTrueVertexParticles = 150;
constexpr int kMaxPFParticles = 100;

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
  void GetEventInfo(art::Event const & evt);
  void GetTruthInfo(art::Event const & evt);
  void FillVertexInfo(art::Event const & evt);
  void FillPFParticleInfo(art::Event const & evt);
  void FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);
  void FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);
  void FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);
  void RunTrackSelection(art::Event const & evt);
  void RunShowerSelection(art::Event const & evt);
  TVector3 ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector);

  ////////////////////////////////////////
  // Trees
  ////////////////////////////////////////
  TTree *fTree;
  TTree *fPOTTree;
  double fPOT;
  ////////////////////////////////////////
  // Event info
  ////////////////////////////////////////
  // Event identification
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  //Detector 
  double fT0;
  //Neutrino 
  int fNuPdg; //Interaction PDG
  int fBeamPdg; //PDG at point of creation
  int fNC;    // 1=is NC, 0=otherwise
  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fTargetZ; //Atomic number of scattering target
  double fQ2; 
  double fENu; 
  double fW; //X-Sec params
  double fX;
  double fY;
  double fNuMomX; //Neutrino momentum
  double fNuMomY;
  double fNuMomZ;
  double fNuMomT;
  double fNuX; //Interaction positions
  double fNuY;
  double fNuZ;
  double fNuT;
  //Outgoing particle count
  int fNPiP; //Number of pi+
  int fNPim; //Number of pi-
  int fNPi0; //Number of pi0
  int fNPhotons; //Number of photons
  int fNProtons; //Number of protons
  int fNNeutrons; //Number of neutrinos
  int fNOther; //Number of other particles
  int fNVertexParticles; //The total number of particles attached to the GHEP vertex
  bool fVertexParticleIsGHEP[kDefMaxNTrueVertexParticles];
  int fVertexParticlePDG[kDefMaxNTrueVertexParticles];
  int fVertexParticleStatus[kDefMaxNTrueVertexParticles];
  int fVertexParticleNChildren[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomX[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomY[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomZ[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomT[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndX[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndY[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndZ[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndT[kDefMaxNTrueVertexParticles];
  //Outgoing Lepton 
  int fLepPDG;
  double fMomLepX;
  double fMomLepY;
  double fMomLepZ;
  double fMomLepT;
  double fLepEndX;
  double fLepEndY;
  double fLepEndZ;
  double fLepEndT;
  double fLepNuAngle;
  //Transverse Mom
  double  fNuMomTranMag;
  double  fTargNuclMomTranMag;
  double  fInitalMomTranMag;
  double  fLepMomTranMag;
  double  fNuclRemMomTranMag;
  double  fFinalMomTranMagNoLepNoRem;
  double  fFinalMomTranMagNoLepWithRem;
  double  fFinalMomTranMagWithLepNoRem;
  double  fFinalMomTranMagWithLepWithRem;
  ////////////////////////////////////////
  // Event-level reco
  ////////////////////////////////////////
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
  int fRecoNuVtxNShowers;
  int fRecoNuVtxNTracks;
  int fRecoNuVtxNChildren;
  double fRecoEventCharge;
  double fNumuRecoMomLep;
  double fNumuRecoEHad;
  double fNumuRecoENu;
  double fNueRecoMomLep;
  double fNueRecoEHad;
  double fNueRecoENu;
  ////////////////////////////////////////
  // PFParticle Info
  ////////////////////////////////////////
  // Truth
  int fRecoPFPTruePDG[kMaxPFParticles];
  bool fRecoPFPTruePrimary[kMaxPFParticles];
  double fRecoPFPTrueMomX[kMaxPFParticles];
  double fRecoPFPTrueMomY[kMaxPFParticles];
  double fRecoPFPTrueMomZ[kMaxPFParticles];
  double fRecoPFPTrueMomT[kMaxPFParticles];
  double fRecoPFPTrueStartX[kMaxPFParticles];
  double fRecoPFPTrueStartY[kMaxPFParticles];
  double fRecoPFPTrueStartZ[kMaxPFParticles];
  double fRecoPFPTrueEndX[kMaxPFParticles];
  double fRecoPFPTrueEndY[kMaxPFParticles];
  double fRecoPFPTrueEndZ[kMaxPFParticles];
  // Reco
  int fNRecoPFPs;
  int fRecoPFPSelf[kMaxPFParticles];
  bool fRecoPFPIsPrimary[kMaxPFParticles];
  int fRecoPFPRecoNHits[kMaxPFParticles];
  double fRecoPFPRecoCompleteness[kMaxPFParticles];
  double fRecoPFPRecoHitPurity[kMaxPFParticles];
  double fRecoPFPRecoCharge[kMaxPFParticles];
  double fRecoPFPRecoVertexX[kMaxPFParticles];
  double fRecoPFPRecoVertexY[kMaxPFParticles];
  double fRecoPFPRecoVertexZ[kMaxPFParticles];
  int fRecoPFPRecoNChildPFP[kMaxPFParticles];
  int fRecoPFPRecoNChildTrackPFP[kMaxPFParticles];
  int fRecoPFPRecoNChildShowerPFP[kMaxPFParticles];
  // DeepPan PID
  double fRecoPFPDeepPanMuVar[kMaxPFParticles];
  double fRecoPFPDeepPanPiVar[kMaxPFParticles];
  double fRecoPFPDeepPanProtonVar[kMaxPFParticles];
  // Ivysaurus PID
  double fRecoPFPIvysaurusMuon[kMaxPFParticles];
  double fRecoPFPIvysaurusProton[kMaxPFParticles];
  double fRecoPFPIvysaurusPion[kMaxPFParticles];
  double fRecoPFPIvysaurusElectron[kMaxPFParticles];
  double fRecoPFPIvysaurusPhoton[kMaxPFParticles];
  double fRecoPFPIvysaurusOther[kMaxPFParticles];
  int fRecoPFPIvysaurusParticleType[kMaxPFParticles];
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
  double fRecoTrackRecoUpstreamX[kMaxPFParticles];
  double fRecoTrackRecoUpstreamY[kMaxPFParticles];
  double fRecoTrackRecoUpstreamZ[kMaxPFParticles];
  double fRecoTrackRecoDownstreamX[kMaxPFParticles];
  double fRecoTrackRecoDownstreamY[kMaxPFParticles];
  double fRecoTrackRecoDownstreamZ[kMaxPFParticles];
  double fRecoTrackRecoEndClosestToVertexX[kMaxPFParticles];
  double fRecoTrackRecoEndClosestToVertexY[kMaxPFParticles];
  double fRecoTrackRecoEndClosestToVertexZ[kMaxPFParticles];
  double fRecoTrackRecoLength[kMaxPFParticles];
  int   fRecoTrackRecoContained[kMaxPFParticles];
  int   fRecoTrackRecoMomMethod[kMaxPFParticles];
  double fRecoTrackRecoMomMCS[kMaxPFParticles];
  double fRecoTrackRecoMomContained[kMaxPFParticles];
  // Pandizzle PID
  float fRecoTrackMichelNHits[kMaxPFParticles];
  float fRecoTrackMichelElectronMVA[kMaxPFParticles];
  float fRecoTrackMichelRecoEnergyPlane2[kMaxPFParticles];
  float fRecoTrackDeflecAngleSD[kMaxPFParticles];
  float fRecoTrackLength[kMaxPFParticles];
  float fRecoTrackEvalRatio[kMaxPFParticles];
  float fRecoTrackConcentration[kMaxPFParticles];
  float fRecoTrackCoreHaloRatio[kMaxPFParticles];
  float fRecoTrackConicalness[kMaxPFParticles];
  float fRecoTrackdEdxStart[kMaxPFParticles];
  float fRecoTrackdEdxEnd[kMaxPFParticles];
  float fRecoTrackdEdxEndRatio[kMaxPFParticles];
  double fRecoTrackPandizzleVar[kMaxPFParticles];
  ////////////////////////////////////////
  // Shower Info
  ////////////////////////////////////////
  double fRecoShowerRecoStartX[kMaxPFParticles];
  double fRecoShowerRecoStartY[kMaxPFParticles];
  double fRecoShowerRecoStartZ[kMaxPFParticles];
  double fRecoShowerRecoDirX[kMaxPFParticles];
  double fRecoShowerRecoDirY[kMaxPFParticles];
  double fRecoShowerRecoDirZ[kMaxPFParticles];
  double fRecoShowerRecoLength[kMaxPFParticles];
  double fRecoShowerRecoOpeningAngle[kMaxPFParticles];
  double fRecoShowerRecodEdx[kMaxPFParticles][3];
  int fRecoShowerRecoBestPlane[kMaxPFParticles];
  double fRecoShowerRecoMom[kMaxPFParticles];
  double fRecoShowerRecoEnergy[kMaxPFParticles][3];
  // Pandrizzle
  double fRecoShowerPandrizzleConnectionBDTScore[kMaxPFParticles];
  double fRecoShowerPandrizzlePathwayLengthMin[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxNPostShowerStartHits[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[kMaxPFParticles];
  double fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxFoundHitRatio[kMaxPFParticles];
  double fRecoShowerPandrizzleMaxInitialGapSize[kMaxPFParticles];
  double fRecoShowerPandrizzleMinLargestProjectedGapSize[kMaxPFParticles];
  double fRecoShowerPandrizzleNViewsWithAmbiguousHits[kMaxPFParticles];
  double fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[kMaxPFParticles];
  double fRecoShowerPandrizzleEvalRatio[kMaxPFParticles];
  double fRecoShowerPandrizzleConcentration[kMaxPFParticles];
  double fRecoShowerPandrizzleCoreHaloRatio[kMaxPFParticles];
  double fRecoShowerPandrizzleConicalness[kMaxPFParticles];
  double fRecoShowerPandrizzledEdxBestPlane[kMaxPFParticles];
  double fRecoShowerPandrizzleDisplacement[kMaxPFParticles];
  double fRecoShowerPandrizzleDCA[kMaxPFParticles];
  double fRecoShowerPandrizzleWideness[kMaxPFParticles];
  double fRecoShowerPandrizzleEnergyDensity[kMaxPFParticles];
  double fRecoShowerPandrizzleBDTMethod[kMaxPFParticles];
  double fRecoShowerEnhancedPandrizzleScore[kMaxPFParticles];
  double fRecoShowerBackupPandrizzleScore[kMaxPFParticles];
  bool   fRecoShowerPandrizzleIsFilled[kMaxPFParticles];
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
  std::string fPOTModuleLabel;
  std::string fCVNModuleLabel;
  ////////////////////////////////////////
  //Algs
  ////////////////////////////////////////
  PandizzleAlg fPandizzleAlg;
  PandrizzleAlg fPandrizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  ctp::CTPHelper fConvTrackPID;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
  ////////////////////////////////////////
  //Tools
  ////////////////////////////////////////
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;
  ////////////////////////////////////////
  // DUNE CVN Scores
  ////////////////////////////////////////
  double fCVNResultNue;
  double fCVNResultNumu;
  double fCVNResultNutau;
  double fCVNResultNC;
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
  fPOTModuleLabel(pset.get<std::string>("POTModuleLabel")),
  fCVNModuleLabel(pset.get<std::string>("CVNModuleLabel")),
  fPandizzleAlg(pset),
  fPandrizzleAlg(pset),
  fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
  fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"), fTrackModuleLabel, fShowerModuleLabel,
  fHitsModuleLabel, fWireModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fRecoModuleLabel),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))}
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
    Reset();

    fRun = evt.run();
    fSubRun = evt.subRun();
    fEvent = evt.event();
    fIsMC = !evt.isRealData();

    GetEventInfo(evt);

    if (fIsMC) 
        GetTruthInfo(evt);

    FillVertexInfo(evt);
    FillPFParticleInfo(evt);
    RunTrackSelection(evt);
    RunShowerSelection(evt);

    fTree->Fill();
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::beginJob()
{
    ///////////////////////////
    // POT tree
    ///////////////////////////
    art::ServiceHandle<art::TFileService> tfs;
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("POT",&fPOT);
    fPOTTree->Branch("Run",&fRun);
    fPOTTree->Branch("SubRun",&fSubRun);

    ///////////////////////////
    // CCNuSelection tree
    ///////////////////////////
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");

    ////////////////////////////
    // Event info
    ////////////////////////////
    fTree->Branch("Run",&fRun);
    fTree->Branch("SubRun",&fSubRun);
    fTree->Branch("Event",&fEvent);

    ////////////////////////////
    // True info
    ////////////////////////////
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
    fTree->Branch("TargetZ",&fTargetZ);
    fTree->Branch("Enu",&fENu);
    fTree->Branch("NuX",&fNuX);
    fTree->Branch("NuY",&fNuY);
    fTree->Branch("NuZ",&fNuZ);
    fTree->Branch("IsMC",&fIsMC);
    fTree->Branch("T0",&fT0);
    fTree->Branch("Q2",&fQ2);
    fTree->Branch("W",&fW);
    fTree->Branch("X",&fX);
    fTree->Branch("Y",&fY);
    fTree->Branch("NuMomX",&fNuMomX);
    fTree->Branch("NuMomY",&fNuMomY);
    fTree->Branch("NuMomZ",&fNuMomZ);
    fTree->Branch("NuMomT",&fNuMomT);
    fTree->Branch("NuT",&fNuT);
    fTree->Branch("NPiP",&fNPiP);
    fTree->Branch("NPim",&fNPim);
    fTree->Branch("NPi0",&fNPi0);
    fTree->Branch("NPhotons",&fNPhotons);
    fTree->Branch("NProton",&fNProtons);
    fTree->Branch("NNeutrons",&fNNeutrons);
    fTree->Branch("NOther",&fNOther);
    fTree->Branch("NVertexParticles",&fNVertexParticles);
    fTree->Branch("VertexParticleIsGHEP",fVertexParticleIsGHEP,"VertexParticleIsGHEP[NVertexParticles]/O");
    fTree->Branch("VertexParticlePDG",fVertexParticlePDG,"VertexParticlePDG[NVertexParticles]/I");
    fTree->Branch("VertexParticleStatus",fVertexParticleStatus,"VertexParticleStatus[NVertexParticles]/I");
    fTree->Branch("VertexParticleNChildren",fVertexParticleNChildren,"VertexParticleNChildren[NVertexParticles]/I");
    fTree->Branch("VertexParticleMomX",fVertexParticleMomX,"VertexParticleMomX[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomY",fVertexParticleMomY,"VertexParticleMomY[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomZ",fVertexParticleMomZ,"VertexParticleMomZ[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomT",fVertexParticleMomT,"VertexParticleMomT[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndX",fVertexParticleEndX,"VertexParticleEndX[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndY",fVertexParticleEndY,"VertexParticleEndY[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndZ",fVertexParticleEndZ,"VertexParticleEndZ[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndT",fVertexParticleEndT,"VertexParticleEndT[NVertexParticles]/D");
    fTree->Branch("LepPDG",&fLepPDG);
    fTree->Branch("MomLepX",&fMomLepX);
    fTree->Branch("MomLepY",&fMomLepY);
    fTree->Branch("MomLepZ",&fMomLepZ);
    fTree->Branch("MomLepT",&fMomLepT);
    fTree->Branch("LepEndX",&fLepEndX);
    fTree->Branch("LepEndY",&fLepEndY);
    fTree->Branch("LepEndZ",&fLepEndZ);
    fTree->Branch("LepEndT",&fLepEndT);
    fTree->Branch("LepNuAngle",&fLepNuAngle);
    fTree->Branch("NuMomTranMag",&fNuMomTranMag);
    fTree->Branch("TargNuclMomTranMag",&fTargNuclMomTranMag);
    fTree->Branch("InitalMomTranMag",&fInitalMomTranMag);
    fTree->Branch("LepMomTranMag",&fLepMomTranMag);
    fTree->Branch("NuclRemMomTranMag",&fNuclRemMomTranMag);
    fTree->Branch("FinalMomTranMagNoLepNoRem",&fFinalMomTranMagNoLepNoRem);
    fTree->Branch("FinalMomTranMagNoLepWithRem",&fFinalMomTranMagNoLepWithRem);
    fTree->Branch("FinalMomTranMagWithLepNoRem",&fFinalMomTranMagWithLepNoRem);
    fTree->Branch("FinalMomTranMagWithLepWithRem",&fFinalMomTranMagWithLepWithRem);

    ////////////////////////////
    // Event-level reco info 
    ////////////////////////////
    fTree->Branch("RecoNuVtxX",&fRecoNuVtxX);
    fTree->Branch("RecoNuVtxY",&fRecoNuVtxY);
    fTree->Branch("RecoNuVtxZ",&fRecoNuVtxZ);
    fTree->Branch("RecoNuVtxNShowers",&fRecoNuVtxNShowers);
    fTree->Branch("RecoNuVtxNTracks",&fRecoNuVtxNTracks);
    fTree->Branch("RecoNuVtxNChildren",&fRecoNuVtxNChildren);
    fTree->Branch("RecoEventCharge",&fRecoEventCharge);
    fTree->Branch("NumuRecoENu",&fNumuRecoENu);
    fTree->Branch("NumuRecoMomLep",&fNumuRecoMomLep);
    fTree->Branch("NumuRecoEHad",&fNumuRecoEHad);
    fTree->Branch("NueRecoENu",&fNueRecoENu);
    fTree->Branch("NueRecoMomLep",&fNueRecoMomLep);
    fTree->Branch("NueRecoEHad",&fNueRecoEHad);

    ////////////////////////////
    // CVN info
    ////////////////////////////
    fTree->Branch("CVNResultNue", &fCVNResultNue);
    fTree->Branch("CVNResultNumu", &fCVNResultNumu);
    fTree->Branch("CVNResultNutau", &fCVNResultNutau);
    fTree->Branch("CVNResultNC", &fCVNResultNC);

    ////////////////////////////
    // PFParticle info
    ////////////////////////////
    fTree->Branch("NRecoPFPs", &fNRecoPFPs);
    // True
    fTree->Branch("RecoPFPTruePDG", fRecoPFPTruePDG, "RecoPFPTruePDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTruePrimary",fRecoPFPTruePrimary,"RecoPFPTruePrimary[NRecoPFPs]/O");
    fTree->Branch("RecoPFPTrueMomX",fRecoPFPTrueMomX,"RecoPFPTrueMomX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomY",fRecoPFPTrueMomY,"RecoPFPTrueMomY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomZ",fRecoPFPTrueMomZ,"RecoPFPTrueMomZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomT",fRecoPFPTrueMomT,"RecoPFPTrueMomT[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartX",fRecoPFPTrueStartX,"RecoPFPTrueStartX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartY",fRecoPFPTrueStartY,"RecoPFPTrueStartY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartZ",fRecoPFPTrueStartZ,"RecoPFPTrueStartZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndX",fRecoPFPTrueEndX,"RecoPFPTrueEndX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndY",fRecoPFPTrueEndY,"RecoPFPTrueEndY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndZ",fRecoPFPTrueEndZ,"RecoPFPTrueEndZ[NRecoPFPs]/D");
    // Reco
    fTree->Branch("RecoPFPSelf", fRecoPFPSelf, "RecoPFPSelf[NRecoPFPs]/I");
    fTree->Branch("RecoPFPIsPrimary", fRecoPFPIsPrimary, "RecoPFPIsPrimary[NRecoPFPs]/O");
    fTree->Branch("RecoPFPRecoNHits",fRecoPFPRecoNHits,"RecoPFPRecoNHits[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoCompleteness",fRecoPFPRecoCompleteness,"RecoPFPRecoCompleteness[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoHitPurity",fRecoPFPRecoHitPurity,"RecoPFPRecoHitPurity[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoCharge",fRecoPFPRecoCharge,"RecoPFPRecoCharge[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexX",fRecoPFPRecoVertexX,"RecoPFPRecoVertexX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexY",fRecoPFPRecoVertexY,"RecoPFPRecoVertexY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexZ",fRecoPFPRecoVertexZ,"RecoPFPRecoVertexZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoNChildPFP",fRecoPFPRecoNChildPFP,"RecoPFPRecoNChildPFP[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoNChildTrackPFP",fRecoPFPRecoNChildTrackPFP,"RecoPFPRecoNChildTrackPFP[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoNChildShowerPFP",fRecoPFPRecoNChildShowerPFP,"RecoPFPRecoNChildShowerPFP[NRecoPFPs]/I");
    // DeepPan
    fTree->Branch("RecoPFPDeepPanMuVar", fRecoPFPDeepPanMuVar, "RecoPFPDeepPanMuVar[NRecoPFPs]/D");
    fTree->Branch("RecoPFPDeepPanPiVar", fRecoPFPDeepPanPiVar, "RecoPFPDeepPanPiVar[NRecoPFPs]/D");
    fTree->Branch("RecoPFPDeepPanProtonVar", fRecoPFPDeepPanProtonVar, "RecoPFPDeepPanProtonVar[NRecoPFPs]/D");
    // Ivysaurus
    fTree->Branch("RecoPFPIvysaurusMuon", fRecoPFPIvysaurusMuon, "RecoPFPIvysaurusMuon[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusProton", fRecoPFPIvysaurusProton, "RecoPFPIvysaurusProton[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusPion", fRecoPFPIvysaurusPion, "RecoPFPIvysaurusPion[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusElectron", fRecoPFPIvysaurusElectron, "RecoPFPIvysaurusElectron[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusPhoton", fRecoPFPIvysaurusPhoton, "RecoPFPIvysaurusPhoton[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusOther", fRecoPFPIvysaurusOther, "RecoPFPIvysaurusOther[NRecoPFPs]/D");
    fTree->Branch("RecoPFPIvysaurusParticleType", fRecoPFPIvysaurusParticleType, "RecoPFPIvysaurusParticleType[NRecoPFPs]/D");

    ////////////////////////////
    // Track Info
    ////////////////////////////
    // Reco
    fTree->Branch("RecoTrackRecoStartX",fRecoTrackRecoStartX,"RecoTrackRecoStartX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartY",fRecoTrackRecoStartY,"RecoTrackRecoStartY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartZ",fRecoTrackRecoStartZ,"RecoTrackRecoStartZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndX",fRecoTrackRecoEndX,"RecoTrackRecoEndX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndY",fRecoTrackRecoEndY,"RecoTrackRecoEndY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndZ",fRecoTrackRecoEndZ,"RecoTrackRecoEndZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamX",fRecoTrackRecoUpstreamX,"RecoTrackRecoUpstreamX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamY",fRecoTrackRecoUpstreamY,"RecoTrackRecoUpstreamY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamZ",fRecoTrackRecoUpstreamZ,"RecoTrackRecoUpstreamZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamX",fRecoTrackRecoDownstreamX,"RecoTrackRecoDownstreamX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamY",fRecoTrackRecoDownstreamY,"RecoTrackRecoDownstreamY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamZ",fRecoTrackRecoDownstreamZ,"RecoTrackRecoDownstreamZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexX",fRecoTrackRecoEndClosestToVertexX,"RecoTrackRecoEndClosestToVertexX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexY",fRecoTrackRecoEndClosestToVertexY,"RecoTrackRecoEndClosestToVertexY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexZ",fRecoTrackRecoEndClosestToVertexZ,"RecoTrackRecoEndClosestToVertexZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoLength",fRecoTrackRecoLength,"RecoTrackRecoLength[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoContained",fRecoTrackRecoContained,"RecoTrackRecoContained[NRecoTracks]/I");
    fTree->Branch("RecoTrackRecoMomMethod",fRecoTrackRecoMomMethod,"RecoTrackRecoMomMethod[NRecoTracks]/I");
    fTree->Branch("RecoTrackRecoMomMCS",fRecoTrackRecoMomMCS,"RecoTrackRecoMomMCS[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomContained",fRecoTrackRecoMomContained,"RecoTrackRecoMomContained[NRecoTracks]/D");
    // Pandizzle
    fTree->Branch("RecoTrackMichelNHits", fRecoTrackMichelNHits, "RecoTrackMichelNHits/D");
    fTree->Branch("RecoTrackMichelElectronMVA", fRecoTrackMichelElectronMVA, "RecoTrackMichelElectronMVA/D");
    fTree->Branch("RecoTrackMichelRecoEnergyPlane2", fRecoTrackMichelRecoEnergyPlane2, "RecoTrackMichelRecoEnergyPlane2/D");
    fTree->Branch("RecoTrackDeflecAngleSD", fRecoTrackDeflecAngleSD, "RecoTrackDeflecAngleSD/D");
    fTree->Branch("RecoTrackLength", fRecoTrackLength, "RecoTrackLength/D");
    fTree->Branch("RecoTrackEvalRatio", fRecoTrackEvalRatio, "RecoTrackEvalRatio/D");
    fTree->Branch("RecoTrackConcentration", fRecoTrackConcentration, "RecoTrackConcentration/D");
    fTree->Branch("RecoTrackCoreHaloRatio", fRecoTrackCoreHaloRatio, "RecoTrackCoreHaloRatio/D");
    fTree->Branch("RecoTrackConicalness", fRecoTrackConicalness, "RecoTrackConicalness/D");
    fTree->Branch("RecoTrackdEdxStart", fRecoTrackdEdxStart, "RecoTrackdEdxStart/D");
    fTree->Branch("RecoTrackdEdxEnd", fRecoTrackdEdxEnd, "RecoTrackdEdxEnd/D");
    fTree->Branch("RecoTrackdEdxEndRatio", fRecoTrackdEdxEndRatio, "RecoTrackdEdxEndRatio/D");

    ///////////////////////////
    // Shower Info
    ///////////////////////////
    // Reco
    fTree->Branch("RecoShowerRecoStartX",fRecoShowerRecoStartX,"RecoShowerRecoStartX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartY",fRecoShowerRecoStartY,"RecoShowerRecoStartY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartZ",fRecoShowerRecoStartZ,"RecoShowerRecoStartZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDirX",fRecoShowerRecoDirX,"RecoShowerRecoDirX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDirY",fRecoShowerRecoDirY,"RecoShowerRecoDirY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDirZ",fRecoShowerRecoDirZ,"RecoShowerRecoDirZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoLength",&fRecoShowerRecoLength,"RecoShowerRecoLength[3]/D");
    fTree->Branch("RecoShowerRecoOpeningAngle",&fRecoShowerRecoOpeningAngle,"RecoShowerRecoOpeningAngle[3]/D");
    fTree->Branch("RecoShowerRecodEdx",fRecoShowerRecodEdx,"RecoShowerRecodEdx[NRecoShowers][3]/D");
    fTree->Branch("RecoShowerRecoBestPlane",&fRecoShowerRecoBestPlane,"RecoShowerRecoBestPlane[3]/I");
    fTree->Branch("RecoShowerRecoMom",fRecoShowerRecoMom,"RecoShowerRecoMom[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEnergy",fRecoShowerRecoEnergy,"RecoShowerRecoEnergy[NRecoShowers][3]/D");
    // Pandrizzle
    fTree->Branch("RecoShowerPandrizzleConnectionBDTScore", &fRecoShowerPandrizzleConnectionBDTScore, "RecoShowerPandrizzleConnectionBDTScore[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzlePathwayLengthMin", &fRecoShowerPandrizzlePathwayLengthMin, "RecoShowerPandrizzlePathwayLengthMin[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D", &fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D, 
                  "RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxNPostShowerStartHits", &fRecoShowerPandrizzleMaxNPostShowerStartHits, "RecoShowerPandrizzleMaxNPostShowerStartHits[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartScatterAngle", &fRecoShowerPandrizzleMaxPostShowerStartScatterAngle, "RecoShowerPandrizzleMaxPostShowerStartScatterAngle[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry, 
                  "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry, 
                  "RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance, 
                  "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius", &fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius, 
                  "RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartOpeningAngle", &fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle, "RecoShowerPandrizzleMaxPostShowerStartOpeningAngle[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxFoundHitRatio", &fRecoShowerPandrizzleMaxFoundHitRatio, "RecoShowerPandrizzleMaxFoundHitRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxInitialGapSize", &fRecoShowerPandrizzleMaxInitialGapSize, "RecoShowerPandrizzleMaxInitialGapSize[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMinLargestProjectedGapSize", &fRecoShowerPandrizzleMinLargestProjectedGapSize, "RecoShowerPandrizzleMinLargestProjectedGapSize[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleNViewsWithAmbiguousHits", &fRecoShowerPandrizzleNViewsWithAmbiguousHits, "RecoShowerPandrizzleNViewsWithAmbiguousHits[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy", &fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy, "RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleEvalRatio",&fRecoShowerPandrizzleEvalRatio,"RecoShowerPandrizzleEvalRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleConcentration",&fRecoShowerPandrizzleConcentration,"RecoShowerPandrizzleConcentration[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleCoreHaloRatio",&fRecoShowerPandrizzleCoreHaloRatio, "RecoShowerPandrizzleCoreHaloRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleConicalness",&fRecoShowerPandrizzleConicalness, "RecoShowerPandrizzleConicalness[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzledEdxBestPlane",&fRecoShowerPandrizzledEdxBestPlane, "RecoShowerPandrizzledEdxBestPlane[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleDisplacement",&fRecoShowerPandrizzleDisplacement, "RecoShowerPandrizzleDisplacement[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleDCA",&fRecoShowerPandrizzleDCA, "RecoShowerPandrizzleDCA[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleWideness",&fRecoShowerPandrizzleWideness, "RecoShowerPandrizzleWideness[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleEnergyDensity",&fRecoShowerPandrizzleEnergyDensity, "RecoShowerPandrizzleEnergyDensity[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleBDTMethod", &fRecoShowerPandrizzleBDTMethod, "RecoShowerPandrizzleBDTMethod[NRecoShowers]/D");
    fTree->Branch("RecoShowerEnhancedPandrizzleScore",&fRecoShowerEnhancedPandrizzleScore, "RecoShowerEnhancedPandrizzleScore[NRecoShowers]/D"); 
    fTree->Branch("RecoShowerBackupPandrizzleScore",&fRecoShowerBackupPandrizzleScore, "RecoShowerBackupPandrizzleScore[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleIsFilled",&fRecoShowerPandrizzleIsFilled, "RecoShowerPandrizzleIsFilled[NRecoShowers]/O");
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::beginSubRun(art::SubRun const & sr)
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::endSubRun(const art::SubRun& sr)
{
  fRun = sr.run();
  fSubRun = sr.subRun();

  // Get POT
  art::Handle<sumdata::POTSummary> potListHandle;

  if (sr.getByLabel(fPOTModuleLabel,potListHandle))
    fPOT = potListHandle->totpot;
  else
    fPOT = 0.;

  if (fPOTTree) fPOTTree->Fill();
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::endJob()
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::Reset()
{
    ////////////////////////////
    // Event info
    ////////////////////////////
    fRun = kDefInt;
    fSubRun = kDefInt;
    fEvent = kDefInt;
    ////////////////////////////
    // True info
    ////////////////////////////
    fIsMC = kDefInt;
    fT0 = kDefDoub;
    fNuPdg = kDefInt; 
    fBeamPdg = kDefInt; 
    fNC = kDefInt;    
    fMode = kDefInt; 
    fTargetZ = kDefInt;
    fQ2 = kDefDoub; 
    fENu = kDefDoub; 
    fW = kDefDoub; 
    fX = kDefDoub;
    fY = kDefDoub;
    fNuMomX = kDefDoub; 
    fNuMomY = kDefDoub;
    fNuMomZ = kDefDoub;
    fNuMomT = kDefDoub;
    fNuX = kDefDoub; 
    fNuY = kDefDoub;
    fNuZ = kDefDoub;
    fNuT = kDefDoub;
    fNPiP = 0;
    fNPim = 0;
    fNPi0 = 0;
    fNPhotons = 0;
    fNProtons = 0;
    fNNeutrons = 0;
    fNOther = 0;
    for (int i = 0; i < kDefMaxNTrueVertexParticles; i++)
    {
        fVertexParticleIsGHEP[i] = 0;
        fVertexParticlePDG[i] = kDefInt;
        fVertexParticleStatus[i] = kDefInt;
        fVertexParticleNChildren[i] = kDefInt;
        fVertexParticleMomX[i] = kDefDoub;
        fVertexParticleMomY[i] = kDefDoub;
        fVertexParticleMomZ[i] = kDefDoub;
        fVertexParticleMomT[i] = kDefDoub;
        fVertexParticleEndX[i] = kDefDoub;
        fVertexParticleEndY[i] = kDefDoub;
        fVertexParticleEndZ[i] = kDefDoub;
        fVertexParticleEndT[i] = kDefDoub;
    }
    fNVertexParticles = 0;
    fLepPDG = kDefInt;
    fMomLepX = kDefDoub;
    fMomLepY = kDefDoub;
    fMomLepZ = kDefDoub;
    fMomLepT = kDefDoub;
    fLepEndX = kDefDoub;
    fLepEndY = kDefDoub;
    fLepEndZ = kDefDoub;
    fLepEndT = kDefDoub;
    fLepNuAngle = kDefDoub;
    fNuMomTranMag = kDefDoub;
    fTargNuclMomTranMag = kDefDoub;
    fInitalMomTranMag = kDefDoub;
    fLepMomTranMag = kDefDoub;
    fNuclRemMomTranMag = kDefDoub;
    fFinalMomTranMagNoLepNoRem = kDefDoub;
    fFinalMomTranMagNoLepWithRem = kDefDoub;
    fFinalMomTranMagWithLepNoRem = kDefDoub;
    fFinalMomTranMagWithLepWithRem = kDefDoub;
    ////////////////////////////
    // Event-level reco info 
    ////////////////////////////
    fRecoNuVtxX = kDefDoub;
    fRecoNuVtxY = kDefDoub;
    fRecoNuVtxZ = kDefDoub;
    fRecoNuVtxNShowers = 0;
    fRecoNuVtxNTracks = 0;
    fRecoNuVtxNChildren = 0;
    fRecoEventCharge = kDefDoub;
    fNumuRecoENu = kDefDoub;
    fNumuRecoMomLep = kDefDoub;
    fNumuRecoEHad = kDefDoub;
    fNueRecoENu = kDefDoub;
    fNueRecoMomLep = kDefDoub;
    fNueRecoEHad = kDefDoub;
    ////////////////////////////
    // CVN info
    ////////////////////////////
    fCVNResultNue = kDefDoub; 
    fCVNResultNumu = kDefDoub;
    fCVNResultNutau = kDefDoub;
    fCVNResultNC = kDefDoub;
    ////////////////////////////
    // PFParticle info
    ////////////////////////////
    fNRecoPFPs = 0;
    for (int i = 0; i < kMaxPFParticles; i++)
    {
        ////////////////////////////  
        // PFParticle stuff
        ////////////////////////////  
        // Truth
        fRecoPFPTruePDG[i] = kDefInt;
        fRecoPFPTruePrimary[i] = false;
        fRecoPFPTrueMomX[i] = kDefDoub;
        fRecoPFPTrueMomY[i] = kDefDoub;
        fRecoPFPTrueMomZ[i] = kDefDoub;
        fRecoPFPTrueMomT[i] = kDefDoub;
        fRecoPFPTrueStartX[i] = kDefDoub;
        fRecoPFPTrueStartY[i] = kDefDoub;
        fRecoPFPTrueStartZ[i] = kDefDoub;
        fRecoPFPTrueEndX[i] = kDefDoub;
        fRecoPFPTrueEndY[i] = kDefDoub;
        fRecoPFPTrueEndZ[i] = kDefDoub;
        // Reco
        fRecoPFPSelf[i] = kDefInt;
        fRecoPFPIsPrimary[i] = false;
        fRecoPFPRecoNHits[i] = kDefInt;
        fRecoPFPRecoCompleteness[i] = kDefDoub;
        fRecoPFPRecoHitPurity[i] = kDefDoub;
        fRecoPFPRecoCharge[i] = kDefDoub;
        fRecoPFPRecoVertexX[i] = kDefDoub;
        fRecoPFPRecoVertexY[i] = kDefDoub;
        fRecoPFPRecoVertexZ[i] = kDefDoub;
        fRecoPFPRecoNChildPFP[i] = kDefInt;
        fRecoPFPRecoNChildTrackPFP[i] = kDefInt;
        fRecoPFPRecoNChildShowerPFP[i] = kDefInt;
        // DeepPan
        fRecoPFPDeepPanMuVar[i] = kDefDoub;
        fRecoPFPDeepPanPiVar[i] = kDefDoub;
        fRecoPFPDeepPanProtonVar[i] = kDefDoub;
        // Ivysaurus
        fRecoPFPIvysaurusMuon[i] = kDefDoub;
        fRecoPFPIvysaurusProton[i] = kDefDoub;
        fRecoPFPIvysaurusPion[i] = kDefDoub;
        fRecoPFPIvysaurusElectron[i] = kDefDoub;
        fRecoPFPIvysaurusPhoton[i] = kDefDoub;
        fRecoPFPIvysaurusOther[i] = kDefDoub;
        fRecoPFPIvysaurusParticleType[i] = kDefInt;
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
        fRecoTrackRecoUpstreamX[i] = kDefDoub;
        fRecoTrackRecoUpstreamY[i] = kDefDoub;
        fRecoTrackRecoUpstreamZ[i] = kDefDoub;
        fRecoTrackRecoDownstreamX[i] = kDefDoub;
        fRecoTrackRecoDownstreamY[i] = kDefDoub;
        fRecoTrackRecoDownstreamZ[i] = kDefDoub;
        fRecoTrackRecoEndClosestToVertexX[i] = kDefDoub;
        fRecoTrackRecoEndClosestToVertexY[i] = kDefDoub;
        fRecoTrackRecoEndClosestToVertexZ[i] = kDefDoub;
        fRecoTrackRecoLength[i] = kDefDoub;
        fRecoTrackRecoContained[i] = kDefInt;
        fRecoTrackRecoMomMethod[i] = kDefInt;
        fRecoTrackRecoMomMCS[i] = kDefDoub;
        fRecoTrackRecoMomContained[i] = kDefDoub;
        // Pandizzle
        fRecoTrackPandizzleVar[i] = kDefDoub;
        fRecoTrackMichelNHits[i] = kDefDoub;
        fRecoTrackMichelElectronMVA[i] = kDefDoub;
        fRecoTrackMichelRecoEnergyPlane2[i] = kDefDoub;
        fRecoTrackDeflecAngleSD[i] = kDefDoub;
        fRecoTrackLength[i] = kDefDoub;
        fRecoTrackEvalRatio[i] = kDefDoub;
        fRecoTrackConcentration[i] = kDefDoub;
        fRecoTrackCoreHaloRatio[i] = kDefDoub;
        fRecoTrackConicalness[i] = kDefDoub;
        fRecoTrackdEdxStart[i] = kDefDoub;
        fRecoTrackdEdxEnd[i] = kDefDoub;
        fRecoTrackdEdxEndRatio[i] = kDefDoub;
        ////////////////////////////
        // Shower stuff
        ////////////////////////////
        // Reco
        fRecoShowerRecoStartX[i] = kDefDoub;
        fRecoShowerRecoStartY[i] = kDefDoub;
        fRecoShowerRecoStartZ[i] = kDefDoub;
        fRecoShowerRecoMom[i] = kDefDoub;
        fRecoShowerRecoDirX[i] = kDefDoub;
        fRecoShowerRecoDirY[i] = kDefDoub;
        fRecoShowerRecoDirZ[i] = kDefDoub;
        fRecoShowerRecoBestPlane[i] = kDefInt;
        fRecoShowerRecoLength[i] = kDefDoub;
        fRecoShowerRecoOpeningAngle[i] = kDefDoub;
        for (int j = 0; j < 3; j++)
        {
            fRecoShowerRecodEdx[i][j] = kDefDoub;
            fRecoShowerRecoEnergy[i][j] = kDefDoub;
        }
        // Pandrizzle
        fRecoShowerPandrizzleConnectionBDTScore[i] = kDefDoub;
        fRecoShowerPandrizzlePathwayLengthMin[i] = kDefDoub;
        fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[i] = kDefDoub;
        fRecoShowerPandrizzleMaxNPostShowerStartHits[i] = kDefDoub;
        fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[i] = kDefDoub;
        fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[i] = kDefDoub;
        fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[i] = kDefDoub;
        fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[i] = kDefDoub;
        fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[i] = kDefDoub;
        fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[i] = kDefDoub;
        fRecoShowerPandrizzleMaxFoundHitRatio[i] = kDefDoub;
        fRecoShowerPandrizzleMaxInitialGapSize[i] = kDefDoub;
        fRecoShowerPandrizzleMinLargestProjectedGapSize[i] = kDefDoub;
        fRecoShowerPandrizzleNViewsWithAmbiguousHits[i] = kDefDoub;
        fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[i] = kDefDoub;
        fRecoShowerPandrizzleEvalRatio[i] = kDefDoub;
        fRecoShowerPandrizzleConcentration[i] = kDefDoub;
        fRecoShowerPandrizzleCoreHaloRatio[i] = kDefDoub;
        fRecoShowerPandrizzleConicalness[i] = kDefDoub;
        fRecoShowerPandrizzledEdxBestPlane[i] = kDefDoub;
        fRecoShowerPandrizzleDisplacement[i] = kDefDoub;
        fRecoShowerPandrizzleDCA[i] = kDefDoub;
        fRecoShowerPandrizzleWideness[i] = kDefDoub;
        fRecoShowerPandrizzleEnergyDensity[i] = kDefDoub;
        fRecoShowerPandrizzleBDTMethod[i] = kDefDoub;
        fRecoShowerEnhancedPandrizzleScore[i] = kDefDoub;
        fRecoShowerBackupPandrizzleScore[i] = kDefDoub;
        fRecoShowerPandrizzleIsFilled[i] = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetEventInfo(art::Event const & evt)
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    // T0
    fT0 = trigger_offset(clockData);

    // Get total event charge
    try
    {
        std::vector<art::Ptr<recob::Hit>> hitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
        fRecoEventCharge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, hitList);
    }
    catch(...)
    {
        return;
    }

    // Get CVN results
    art::Handle<std::vector<cvn::Result>> cvnResult;
    evt.getByLabel(fCVNModuleLabel, cvnResult);

    if (!cvnResult->empty()) 
    {
        fCVNResultNue = (*cvnResult)[0].GetNueProbability();
        fCVNResultNumu = (*cvnResult)[0].GetNumuProbability();
        fCVNResultNutau = (*cvnResult)[0].GetNutauProbability();
        fCVNResultNC = (*cvnResult)[0].GetNCProbability();
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetTruthInfo(art::Event const & evt)
{
  //Get the generator record
  art::Handle<std::vector<simb::MCTruth>> mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle)){
    art::fill_ptr_vector(mcList, mcTruthListHandle);
  }

  if (mcList.size() == 0)
    return;

  if (mcList.size() > 1)
    mf::LogWarning("CCNuSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one.";

  //Get the flux record
  art::Handle<std::vector<simb::MCFlux>> mcFluxListHandle;
  std::vector<art::Ptr<simb::MCFlux>> mcFlux;
  if (evt.getByLabel(fNuGenModuleLabel, mcFluxListHandle)){
    art::fill_ptr_vector(mcFlux, mcFluxListHandle);
  }

  //need the assns for later
  art::FindManyP<simb::MCParticle> fmpt(mcTruthListHandle, evt, fLargeantModuleLabel);

  art::Ptr<simb::MCTruth> mcTruth = mcList.at(0);

  if (mcTruth->Origin() != simb::kBeamNeutrino) 
  {
    mf::LogWarning("CCNuSelection") << "Origin for this event is " << mcTruth->Origin() << " and not simb::kBeamNeutrino (" << simb::kBeamNeutrino<<")";
    return;
  }
 
  // Neutrino
  fNuPdg = mcTruth->GetNeutrino().Nu().PdgCode();

  if (mcFluxListHandle.isValid()) 
    fBeamPdg  = mcFlux[0]->fntype;

  fNC = mcTruth->GetNeutrino().CCNC();
  fMode = mcTruth->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  fTargetZ = mcTruth->GetNeutrino().Target()%100000000/10000;
  fENu = mcTruth->GetNeutrino().Nu().E();
  fQ2 = mcTruth->GetNeutrino().QSqr();
  fW = mcTruth->GetNeutrino().W();
  fX = mcTruth->GetNeutrino().X();
  fY = mcTruth->GetNeutrino().Y();
  fNuX = mcTruth->GetNeutrino().Nu().Vx();
  fNuY = mcTruth->GetNeutrino().Nu().Vy();
  fNuZ = mcTruth->GetNeutrino().Nu().Vz();
  fNuT = mcTruth->GetNeutrino().Nu().T();
  fNuMomX = mcTruth->GetNeutrino().Nu().Momentum().X();
  fNuMomY = mcTruth->GetNeutrino().Nu().Momentum().Y();
  fNuMomZ = mcTruth->GetNeutrino().Nu().Momentum().Z();
  fNuMomT = mcTruth->GetNeutrino().Nu().Momentum().T();

  //Leading lepton
  fLepPDG = mcTruth->GetNeutrino().Lepton().PdgCode();
  fMomLepX = mcTruth->GetNeutrino().Lepton().Momentum().X();
  fMomLepY = mcTruth->GetNeutrino().Lepton().Momentum().Y();
  fMomLepZ = mcTruth->GetNeutrino().Lepton().Momentum().Z();
  fMomLepT = mcTruth->GetNeutrino().Lepton().Momentum().T();
  fLepEndX = mcTruth->GetNeutrino().Lepton().EndPosition().X();
  fLepEndY = mcTruth->GetNeutrino().Lepton().EndPosition().Y();
  fLepEndZ = mcTruth->GetNeutrino().Lepton().EndPosition().Z();
  fLepEndY = mcTruth->GetNeutrino().Lepton().EndPosition().T();
  fLepNuAngle = mcTruth->GetNeutrino().Nu().Momentum().Vect().Angle(mcTruth->GetNeutrino().Lepton().Momentum().Vect());

  //Vertex particles
  const std::vector<art::Ptr<simb::MCParticle>> associated_particles = fmpt.at(mcTruth.key());

  fNVertexParticles = 0;

  for (unsigned int i_part = 0; i_part < associated_particles.size(); i_part++)
  {
    art::Ptr<simb::MCParticle> particle = associated_particles[i_part];
    if (particle->StatusCode() != 1) continue; //count tracked particles
    if (particle->Mother() != 0) continue; //count primary particles
    fVertexParticleIsGHEP[fNVertexParticles] = 0;
    fVertexParticlePDG[fNVertexParticles] = particle->PdgCode();;
    fVertexParticleStatus[fNVertexParticles] = particle->StatusCode();
    fVertexParticleNChildren[fNVertexParticles] = particle->NumberDaughters();
    fVertexParticleMomX[fNVertexParticles] = particle->Momentum(0).X();
    fVertexParticleMomY[fNVertexParticles] = particle->Momentum(0).Y();
    fVertexParticleMomZ[fNVertexParticles] = particle->Momentum(0).Z();
    fVertexParticleMomT[fNVertexParticles] = particle->Momentum(0).T();
    fVertexParticleEndX[fNVertexParticles] = particle->EndPosition().X();
    fVertexParticleEndY[fNVertexParticles] = particle->EndPosition().Y();
    fVertexParticleEndZ[fNVertexParticles] = particle->EndPosition().Z();
    fVertexParticleEndT[fNVertexParticles] = particle->EndPosition().T();
    fNVertexParticles++;
  }

  //Loop over the final state particles from the ghep vertex
  for (int i_part = 0; i_part < mcTruth->NParticles(); i_part++)
  {
    const simb::MCParticle& vertex_particle = mcTruth->GetParticle(i_part);
    int pdg = vertex_particle.PdgCode();
    fVertexParticleIsGHEP[fNVertexParticles] = 1;
    fVertexParticlePDG[fNVertexParticles] = pdg;
    fVertexParticleStatus[fNVertexParticles] = vertex_particle.StatusCode();
    fVertexParticleNChildren[fNVertexParticles] = vertex_particle.NumberDaughters();
    fVertexParticleMomX[fNVertexParticles] = vertex_particle.Momentum(0).X();
    fVertexParticleMomY[fNVertexParticles] = vertex_particle.Momentum(0).Y();
    fVertexParticleMomZ[fNVertexParticles] = vertex_particle.Momentum(0).Z();
    fVertexParticleMomT[fNVertexParticles] = vertex_particle.Momentum(0).T();
    fVertexParticleEndX[fNVertexParticles] = vertex_particle.EndPosition().X();
    fVertexParticleEndY[fNVertexParticles] = vertex_particle.EndPosition().Y();
    fVertexParticleEndZ[fNVertexParticles] = vertex_particle.EndPosition().Z();
    fVertexParticleEndT[fNVertexParticles] = vertex_particle.EndPosition().T();
    fNVertexParticles++;
    if (!(vertex_particle.StatusCode() == 1)) continue;
    if (pdg >= 2000000000) continue;
    if (std::abs(pdg) >= 11 && std::abs(pdg) <= 16) continue;
    if (pdg==211) fNPiP++;
    else if (pdg==-211) fNPim++;
    else if (pdg==111) fNPi0++;
    else if (pdg==22) fNPhotons++;
    else if (pdg==2212) fNProtons++;
    else if (pdg==2112) fNNeutrons++;
    else fNOther++;
  }

  if (fNVertexParticles > kDefMaxNTrueVertexParticles)
  {
    std::cout << "CCNuSelection::GetTruthInfo VERTEX ARRAY IS GOING TO BE OVERFILLED - VERY BAD" << std::endl;
    throw;
  }

  //do some transverse momentum stuff
  TVector3 beam_axis(0, 0, 1);
  beam_axis.RotateX(-0.101);
  TVector3 nu_mom_vect = mcTruth->GetNeutrino().Nu().Momentum().Vect();
  TVector3 nu_mom_tran_vect = ProjectVectorOntoPlane(nu_mom_vect, beam_axis);
  fNuMomTranMag = nu_mom_tran_vect.Mag(); 

  TVector3 total_initial_mom_tran_vect;
  total_initial_mom_tran_vect += nu_mom_tran_vect;

  TVector3 total_final_mom_tran_vect_nolep_norem;
  TVector3 total_final_mom_tran_vect_nolep_withrem;
  TVector3 total_final_mom_tran_vect_withlep_norem;
  TVector3 total_final_mom_tran_vect_withlep_withrem;
  
  //Loop over the particles
  for (int i_part = 0; i_part < mcTruth->NParticles(); i_part++){
    const simb::MCParticle& particle = mcTruth->GetParticle(i_part);
    int status = particle.StatusCode();
    if (status == 11){ //Nucleon target
      TVector3 target_nucleon_mom_vect = particle.Momentum(0).Vect();
      TVector3 target_nucleon_mom_tran_vect = ProjectVectorOntoPlane(target_nucleon_mom_vect, beam_axis); 
      total_initial_mom_tran_vect += target_nucleon_mom_tran_vect;
      fTargNuclMomTranMag = target_nucleon_mom_tran_vect.Mag();
    }
    else if (status == 15){ //final nuclear remnant
      TVector3 nuclear_rem_mom_vect = particle.Momentum(0).Vect();
      TVector3 nuclear_rem_mom_tran_vect = ProjectVectorOntoPlane(nuclear_rem_mom_vect, beam_axis);
      total_final_mom_tran_vect_nolep_withrem += nuclear_rem_mom_tran_vect;
      total_final_mom_tran_vect_withlep_withrem += nuclear_rem_mom_tran_vect;
      fNuclRemMomTranMag = nuclear_rem_mom_tran_vect.Mag();
    }
    else if (status == 1){ //final state particle
      TVector3 fin_state_mom_vect = particle.Momentum(0).Vect();
      TVector3 fin_state_mom_tran_vect = ProjectVectorOntoPlane(fin_state_mom_vect, beam_axis);
      int pdg = particle.PdgCode();
      if (std::abs(pdg) >= 11 && std::abs(pdg) <= 16){
        total_final_mom_tran_vect_withlep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_withlep_norem += fin_state_mom_tran_vect;
        fLepMomTranMag = fin_state_mom_tran_vect.Mag();
      }
      else {
        total_final_mom_tran_vect_withlep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_withlep_norem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_nolep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_nolep_norem += fin_state_mom_tran_vect;
      }
    }
  }

  //Get all of the mom mags now
  fInitalMomTranMag = total_initial_mom_tran_vect.Mag();
  fFinalMomTranMagNoLepNoRem = total_final_mom_tran_vect_nolep_norem.Mag();
  fFinalMomTranMagNoLepWithRem = total_final_mom_tran_vect_nolep_withrem.Mag();
  fFinalMomTranMagWithLepNoRem = total_final_mom_tran_vect_withlep_norem.Mag();
  fFinalMomTranMagWithLepWithRem = total_final_mom_tran_vect_withlep_withrem.Mag();
}

///////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillVertexInfo(art::Event const & evt)
{
  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
    return;

  art::Ptr<recob::PFParticle> nu_pfp = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fRecoModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nu_pfp, evt, fRecoModuleLabel);

  fRecoNuVtxNChildren = nuChildren.size();

  for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
  {
    int pdg = nuChild->PdgCode();

    if (pdg == 11) 
      fRecoNuVtxNShowers++;
    else if (pdg == 13) 
      fRecoNuVtxNTracks++;
  }

  try
  {
      art::Ptr<recob::Vertex> nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nu_pfp, evt, fRecoModuleLabel);

      fRecoNuVtxX = nuVertex->position().X();
      fRecoNuVtxY = nuVertex->position().Y();
      fRecoNuVtxZ = nuVertex->position().Z();
  }
  catch(...)
  {
      return;
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillPFParticleInfo(art::Event const & evt)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
        return;

    art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fRecoModuleLabel);
    std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fRecoModuleLabel); 
    std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);
    std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);

    int pfpCounter = -1;

    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        pfpCounter++;

        fRecoPFPSelf[pfpCounter] = pfp->Self();

        if (pfpCounter == kMaxPFParticles)
            break;

        for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
        {
            if (nuChild->Self() == pfp->Self())
            {
                fRecoPFPIsPrimary[pfpCounter] = true;
                break;
            }
        }

        const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fRecoModuleLabel);
        fRecoPFPRecoNHits[pfpCounter] = pfpHits.size();

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
        fRecoPFPRecoCharge[pfpCounter]  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, pfpHits);

        try
        {
            art::Ptr<recob::Vertex> recoVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, fRecoModuleLabel);

            fRecoPFPRecoVertexX[pfpCounter] = recoVertex->position().X();
            fRecoPFPRecoVertexY[pfpCounter] = recoVertex->position().Y();
            fRecoPFPRecoVertexZ[pfpCounter] = recoVertex->position().Z();
        }
        catch(...)
        {
        }

        // Child particle info
        FillChildPFPInformation(pfp, evt, fRecoPFPRecoNChildPFP[pfpCounter], fRecoPFPRecoNChildTrackPFP[pfpCounter], fRecoPFPRecoNChildShowerPFP[pfpCounter]);

        int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);
        fRecoPFPRecoCompleteness[pfpCounter] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, pfpHits, eventHitList, g4id);
        fRecoPFPRecoHitPurity[pfpCounter] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, pfpHits, g4id);

        if (TruthMatchUtils::Valid(g4id))
        {
            art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
            const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);

            if (matched_mcparticle)
            {
                fRecoPFPTruePDG[pfpCounter] = matched_mcparticle->PdgCode();

                if (matched_mcparticle->Mother() == 0) 
                    fRecoPFPTruePrimary[pfpCounter] = true;
                else 
                    fRecoPFPTruePrimary[pfpCounter] = false;

                fRecoPFPTrueMomX[pfpCounter] = matched_mcparticle->Momentum().X();
                fRecoPFPTrueMomY[pfpCounter] = matched_mcparticle->Momentum().Y();
                fRecoPFPTrueMomZ[pfpCounter] = matched_mcparticle->Momentum().Z();
                fRecoPFPTrueStartX[pfpCounter] = matched_mcparticle->Position(0).X();
                fRecoPFPTrueStartY[pfpCounter] = matched_mcparticle->Position(0).Y();
                fRecoPFPTrueStartZ[pfpCounter] = matched_mcparticle->Position(0).Z();
                fRecoPFPTrueEndX[pfpCounter] = matched_mcparticle->EndPosition().X();
                fRecoPFPTrueEndY[pfpCounter] = matched_mcparticle->EndPosition().Y();
                fRecoPFPTrueEndZ[pfpCounter] = matched_mcparticle->EndPosition().Z();
            }
        }

        // Ivysaurus stuff

        // DeepPan stuff
        ctp::CTPResult deepPanPIDResult = fConvTrackPID.RunConvolutionalTrackPID(pfp, evt);
        if (deepPanPIDResult.IsValid())
        {
            fRecoPFPDeepPanMuVar[pfpCounter] = deepPanPIDResult.GetMuonScore();
            fRecoPFPDeepPanPiVar[pfpCounter] = deepPanPIDResult.GetPionScore();
            fRecoPFPDeepPanProtonVar[pfpCounter] = deepPanPIDResult.GetProtonScore();
        }

        // Fill the track information
        FillRecoTrackInfo(evt, pfp, pfpCounter);
        FillRecoShowerInfo(evt, pfp, pfpCounter);
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, 
    int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp)
{
    n_child_pfp = 0;
    n_child_track_pfp = 0;
    n_child_shower_pfp = 0;

    std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfp, evt, fRecoModuleLabel);

    for (art::Ptr<recob::PFParticle> childPFP : childPFPs)
    {
        int pdg = childPFP->PdgCode();

        if (pdg == 13)
            n_child_track_pfp++;
        else if (pdg == 11)
            n_child_shower_pfp++;
        else 
            std::cout << "FillChildPFPInformation: found a child PFP with an unexpected pdg code: " << pdg << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpCounter)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fRecoModuleLabel, fTrackModuleLabel))
        return;

    art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, fRecoModuleLabel, fTrackModuleLabel);

    // General track variables
    recob::Track::Point_t trackStart, trackEnd;
    std::tie(trackStart, trackEnd) = track->Extent(); 
    fRecoTrackRecoStartX[pfpCounter] = trackStart.X();
    fRecoTrackRecoStartY[pfpCounter] = trackStart.Y();
    fRecoTrackRecoStartZ[pfpCounter] = trackStart.Z();
    fRecoTrackRecoEndX[pfpCounter] = trackEnd.X();
    fRecoTrackRecoEndY[pfpCounter] = trackEnd.Y();
    fRecoTrackRecoEndZ[pfpCounter] = trackEnd.Z();

    if (fRecoTrackRecoEndZ[pfpCounter] > fRecoTrackRecoStartZ[pfpCounter])
    {
        fRecoTrackRecoUpstreamX[pfpCounter] = fRecoTrackRecoStartX[pfpCounter];
        fRecoTrackRecoUpstreamY[pfpCounter] = fRecoTrackRecoStartY[pfpCounter];
        fRecoTrackRecoUpstreamZ[pfpCounter] = fRecoTrackRecoStartZ[pfpCounter];
        fRecoTrackRecoDownstreamX[pfpCounter] = fRecoTrackRecoEndX[pfpCounter];
        fRecoTrackRecoDownstreamY[pfpCounter] = fRecoTrackRecoEndY[pfpCounter];
        fRecoTrackRecoDownstreamZ[pfpCounter] = fRecoTrackRecoEndZ[pfpCounter];
    }
    else
    {
        fRecoTrackRecoDownstreamX[pfpCounter] = fRecoTrackRecoStartX[pfpCounter];
        fRecoTrackRecoDownstreamY[pfpCounter] = fRecoTrackRecoStartY[pfpCounter];
        fRecoTrackRecoDownstreamZ[pfpCounter] = fRecoTrackRecoStartZ[pfpCounter];
        fRecoTrackRecoUpstreamX[pfpCounter] = fRecoTrackRecoEndX[pfpCounter];
        fRecoTrackRecoUpstreamY[pfpCounter] = fRecoTrackRecoEndY[pfpCounter];
        fRecoTrackRecoUpstreamZ[pfpCounter] = fRecoTrackRecoEndZ[pfpCounter];
    }

    fRecoTrackRecoLength[pfpCounter] = track->Length();

    TVector3 upstream_end(fRecoTrackRecoUpstreamX[pfpCounter], fRecoTrackRecoUpstreamY[pfpCounter], fRecoTrackRecoUpstreamZ[pfpCounter]);
    TVector3 downstream_end(fRecoTrackRecoDownstreamX[pfpCounter], fRecoTrackRecoDownstreamY[pfpCounter], fRecoTrackRecoDownstreamZ[pfpCounter]);
    TVector3 vertex_pos(fRecoPFPRecoVertexX[pfpCounter], fRecoPFPRecoVertexY[pfpCounter], fRecoPFPRecoVertexZ[pfpCounter]);

    if ((vertex_pos - upstream_end).Mag() < (vertex_pos - downstream_end).Mag())
    {
        fRecoTrackRecoEndClosestToVertexX[pfpCounter] = fRecoTrackRecoUpstreamX[pfpCounter];
        fRecoTrackRecoEndClosestToVertexY[pfpCounter] = fRecoTrackRecoUpstreamY[pfpCounter];
        fRecoTrackRecoEndClosestToVertexZ[pfpCounter] = fRecoTrackRecoUpstreamZ[pfpCounter];
    }
    else
    {
        fRecoTrackRecoEndClosestToVertexX[pfpCounter] = fRecoTrackRecoDownstreamX[pfpCounter];
        fRecoTrackRecoEndClosestToVertexY[pfpCounter] = fRecoTrackRecoDownstreamY[pfpCounter];
        fRecoTrackRecoEndClosestToVertexZ[pfpCounter] = fRecoTrackRecoDownstreamZ[pfpCounter];
    }

    // Fill momentum variables
    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(track, evt)));
    fRecoTrackRecoContained[pfpCounter] = energyRecoHandle->longestTrackContained;
    fRecoTrackRecoMomMethod[pfpCounter] = energyRecoHandle->trackMomMethod;

    if (energyRecoHandle->trackMomMethod == 1)
        fRecoTrackRecoMomContained[pfpCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    else if (energyRecoHandle->trackMomMethod == 0)
        fRecoTrackRecoMomMCS[pfpCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    // Pandizzle variables
    FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(track, evt));
    fRecoTrackMichelNHits[pfpCounter] = (float)pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits);
    fRecoTrackMichelElectronMVA[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA);
    fRecoTrackMichelRecoEnergyPlane2[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2);
    fRecoTrackDeflecAngleSD[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
    fRecoTrackLength[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
    fRecoTrackEvalRatio[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
    fRecoTrackConcentration[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
    fRecoTrackCoreHaloRatio[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
    fRecoTrackConicalness[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
    fRecoTrackdEdxStart[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
    fRecoTrackdEdxEnd[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
    fRecoTrackdEdxEndRatio[pfpCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);
    fRecoTrackPandizzleVar[pfpCounter] = pandizzleRecord.GetMVAScore();
}

///////////////////////////////////////////////////////////////
// Do different selections here...
void FDSelection::CCNuSelection::RunTrackSelection(art::Event const & evt)
{
  // Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(evt);

  if (!sel_track.isAvailable()) 
  {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no track returned from selection" << std::endl; 
    return;
  }

  // Fill neutrino energy variables
  // Use selected track to get neutrino energy
  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_track, evt)));
  fNumuRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNumuRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fNumuRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpCounter)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel))
        return;

    art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel);

    // General shower variables
    fRecoShowerRecoDirX[pfpCounter] = shower->Direction().X();
    fRecoShowerRecoDirY[pfpCounter] = shower->Direction().Y();
    fRecoShowerRecoDirZ[pfpCounter] = shower->Direction().Z();
    fRecoShowerRecoStartX[pfpCounter] = shower->ShowerStart().X();
    fRecoShowerRecoStartY[pfpCounter] = shower->ShowerStart().Y();
    fRecoShowerRecoStartZ[pfpCounter] = shower->ShowerStart().Z();
    fRecoShowerRecoBestPlane[pfpCounter] = shower->best_plane();
    fRecoShowerRecoLength[pfpCounter] = shower->Length();
    fRecoShowerRecoOpeningAngle[pfpCounter] = shower->OpenAngle();

    // Momentum and energy
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(shower, evt)));
    fRecoShowerRecoMom[pfpCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    if (shower->dEdx().size() > 0)
    {
      for (int i_plane = 0; i_plane < 3; i_plane++)
      {
        fRecoShowerRecodEdx[pfpCounter][i_plane] = shower->dEdx()[i_plane];
        fRecoShowerRecoEnergy[pfpCounter][i_plane] = shower->Energy()[i_plane];
      }
    }

    // Pandrizzle
    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower, evt));

    fRecoShowerPandrizzleEvalRatio[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
    fRecoShowerPandrizzleConcentration[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
    fRecoShowerPandrizzleCoreHaloRatio[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
    fRecoShowerPandrizzleConicalness[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
    fRecoShowerPandrizzledEdxBestPlane[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
    fRecoShowerPandrizzleDisplacement[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
    fRecoShowerPandrizzleDCA[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
    fRecoShowerPandrizzleWideness[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
    fRecoShowerPandrizzleEnergyDensity[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);
    fRecoShowerPandrizzlePathwayLengthMin[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
    fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
    fRecoShowerPandrizzleMaxNPostShowerStartHits[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
    fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[pfpCounter] = 
      pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
    fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
    fRecoShowerPandrizzleMaxFoundHitRatio[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
    fRecoShowerPandrizzleMaxInitialGapSize[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
    fRecoShowerPandrizzleMinLargestProjectedGapSize[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
    fRecoShowerPandrizzleNViewsWithAmbiguousHits[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
    fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);
    fRecoShowerPandrizzleBDTMethod[pfpCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);

    float pandrizzleScore(pandrizzleRecord.GetMVAScore());
    fRecoShowerBackupPandrizzleScore[pfpCounter] = (std::fabs(fRecoShowerPandrizzleBDTMethod[pfpCounter] - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerEnhancedPandrizzleScore[pfpCounter] = (std::fabs(fRecoShowerPandrizzleBDTMethod[pfpCounter] - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerPandrizzleIsFilled[pfpCounter] = pandrizzleRecord.IsFilled();
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunShowerSelection(art::Event const & evt)
{
  // Get the selected shower
  art::Ptr<recob::Shower> sel_shower = fRecoShowerSelector->FindSelectedShower(evt);

  if (!sel_shower.isAvailable()) 
  {
    std::cout << "FDSelection::CCNuSelection::RunShowerSelection - no shower selected by tool" << std::endl;
    return;
  }

  // Momentum and energy
  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_shower, evt)));
  fNueRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNueRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fNueRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
}

//////////////////////////////////////////////////////////////////

TVector3 FDSelection::CCNuSelection::ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector){
  TVector3 projected_vector = vector_to_project - (vector_to_project.Dot(plane_norm_vector) / (plane_norm_vector.Mag() * plane_norm_vector.Mag()))*plane_norm_vector;
  return projected_vector;
}

//////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(FDSelection::CCNuSelection)
