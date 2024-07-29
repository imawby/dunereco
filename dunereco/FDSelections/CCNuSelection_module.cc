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
#include "dunereco/TrackPID/products/CTPResult.h"
#include "dunereco/Ivysaurus/TensorFlow/IvysaurusEvaluator.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

//Custom
#include "dunereco/FDSelections/EnergyReco/TrueEnergyCalc.h"
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "dunereco/FDSelections/HierarchyUtils.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -999;
constexpr double kDefDoub = -999.0;
constexpr int kDefMaxNTrueVertexParticles = 150;
constexpr int kMaxPFParticles = 100;
constexpr int kMaxParentChildLinks = (kMaxPFParticles * kMaxPFParticles);

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
  void GetEventInfo(art::Event const & evt);
  void GetTruthInfo(art::Event const & evt);
  void FillVertexInfo(art::Event const & evt);
  void FillPFParticleInfo(art::Event const & evt);
  void FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);
  void FillHierarchyInfo(art::Event const & evt);
  void SetRecoGenerationInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, int &recoGeneration, int &recoParentSelf, int &recoParentPDG);
  void SetTrueGenerationInfo(art::Event const & evt, const int pfpIndex, const bool visibleMode, int &trueGeneration, int &trueParentTrackID, int &trueParentPDG);
  void FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);
  void FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, const int pfpCounter);
  void FillParentChildLinkInfo(art::Event const & evt);
  void FillTrueParentChildLinkInfo(const int linkIndex, const int parentPFPIndex, const int childPFPIndex);
  void FillRecoParentChildLinkInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> childPFP, art::Ptr<recob::PFParticle> parentPFP, const int linkIndex);
  void RunTrackSelection(art::Event const & evt);
  void RunPandizzleTrackSelection();
  void RunDeepPanTrackSelection();
  void RunIvysaurusTrackSelection();
  void RunLongestLengthTrackSelection();
  void RunShowerSelection(art::Event const & evt);
  void RunPandrizzleShowerSelection();
  void RunIvysaurusShowerSelection();
  void RunHighestEnergyShowerSelection();
  TVector3 ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector);

  ////////////////////////////////////////
  // Pandora maps
  ////////////////////////////////////////
  lar_pandora::MCParticleMap fMCParticleMap;   // Linking TrackID -> MCParticle
  lar_pandora::PFParticleMap fPFPMap;          // Linking Self() -> PFParticle
  std::vector<int> fReconstructedMCParticles;  // TrackIDs of reco'd MCParticles
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
  int fNuTrackID;
  int fNC;    // 1=is NC, 0=otherwise
  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fTargetZ; //Atomic number of scattering target
  double fQ2; 
  double fENu; 
  double fTrueNuEnergyEDep; 
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
  //Outgoing Lepton 
  int fLepPDG;
  double fLepEnergy;
  double fTrueLepEnergyEDep;
  double fMomLepX;
  double fMomLepY;
  double fMomLepZ;
  double fMomLepT;
  double fLepEndX;
  double fLepEndY;
  double fLepEndZ;
  double fLepEndT;
  double fLepNuAngle;
  int fLepHasMichelDecay;
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
  std::vector<float> fTrueEnergyDepX_Inner;
  std::vector<float> fTrueEnergyDepY_Inner;
  std::vector<float> fTrueEnergyDepZ_Inner;
  std::vector<int> fIsEdepInAPA;
  std::vector<float> fAPALowY;
  std::vector<float> fAPAHighY;
  std::vector<float> fAPALowZ;
  std::vector<float> fAPAHighZ;
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
  double fNumuRecoMomLep[kMaxPFParticles];
  double fNumuRecoEHad[kMaxPFParticles];
  double fNumuRecoENu[kMaxPFParticles];
  double fNueRecoMomLep[kMaxPFParticles];
  double fNueRecoEHad[kMaxPFParticles];
  double fNueRecoENu[kMaxPFParticles];
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
  double fRecoPFPTrueEnergyEDep[kMaxPFParticles];
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
  double fRecoPFPTrackShowerScore[kMaxPFParticles];
  int fRecoPFPTrackShowerPDG[kMaxPFParticles];
  int fRecoPFPRecoNHits[kMaxPFParticles];
  std::vector<std::vector<double>> fRecoPFPSpacepointX;
  std::vector<std::vector<double>> fRecoPFPSpacepointY;
  std::vector<std::vector<double>> fRecoPFPSpacepointZ;
  int fRecoPFPRecoGeneration[kMaxPFParticles];
  int fRecoPFPRecoParentSelf[kMaxPFParticles];
  int fRecoPFPRecoParentPDG[kMaxPFParticles];
  int fRecoPFPShowerFitSuccess[kMaxPFParticles];
  int fRecoPFPTrackFitSuccess[kMaxPFParticles];
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
  double fRecoTrackRecoMomRange[kMaxPFParticles];
  // Pandizzle PID
  double fRecoTrackMichelNHits[kMaxPFParticles];
  double fRecoTrackMichelElectronMVA[kMaxPFParticles];
  double fRecoTrackMichelRecoEnergyPlane2[kMaxPFParticles];
  double fRecoTrackDeflecAngleSD[kMaxPFParticles];
  double fRecoTrackLength[kMaxPFParticles];
  double fRecoTrackEvalRatio[kMaxPFParticles];
  double fRecoTrackConcentration[kMaxPFParticles];
  double fRecoTrackCoreHaloRatio[kMaxPFParticles];
  double fRecoTrackConicalness[kMaxPFParticles];
  double fRecoTrackdEdxStart[kMaxPFParticles];
  double fRecoTrackdEdxEnd[kMaxPFParticles];
  double fRecoTrackdEdxEndRatio[kMaxPFParticles];
  double fRecoTrackPandizzleVar[kMaxPFParticles];
  // Selected info
  int fSelTrackPandizzleSelf; 
  int fSelTrackPandizzleIndex;
  int fSelTrackDeepPanSelf; 
  int fSelTrackDeepPanIndex;
  int fSelTrackIvysaurusSelf; 
  int fSelTrackIvysaurusIndex;
  int fSelTrackLongestLengthSelf; 
  int fSelTrackLongestLengthIndex;
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
  double fRecoShowerPandrizzleModularPathwayLength[kMaxPFParticles];
  double fRecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance[kMaxPFParticles];
  double fRecoShowerPandrizzleModularMaxNShowerHits[kMaxPFParticles];
  double fRecoShowerPandrizzleBDTMethod[kMaxPFParticles];
  double fRecoShowerEnhancedPandrizzleScore[kMaxPFParticles];
  double fRecoShowerBackupPandrizzleScore[kMaxPFParticles];
  bool   fRecoShowerPandrizzleIsFilled[kMaxPFParticles];
  // Selected info
  int fSelShowerPandrizzleSelf; 
  int fSelShowerPandrizzleIndex;
  int fSelShowerIvysaurusSelf; 
  int fSelShowerIvysaurusIndex;
  int fSelShowerHighestEnergySelf; 
  int fSelShowerHighestEnergyIndex;
  ////////////////////////////////////////
  // Hierarchy Info
  ////////////////////////////////////////
  // Truth
  bool fTrueParentChildLink[kMaxParentChildLinks];
  // Reco
  int fNParentChildLinks;
  // Parent information
  double fParentTrackScore[kMaxParentChildLinks];
  double fParentBraggVariable[kMaxParentChildLinks];
  double fParentEndRegionNHits[kMaxParentChildLinks];
  double fParentEndRegionNParticles[kMaxParentChildLinks];
  double fParentEndRegionRToWall[kMaxParentChildLinks];
  // Edge information
  double fParentPFPIndex[kMaxParentChildLinks];
  double fChildPFPIndex[kMaxParentChildLinks];
  double fVertexSeparation[kMaxParentChildLinks];
  double fSeparationU[kMaxParentChildLinks];
  double fSeparationV[kMaxParentChildLinks];
  double fSeparationW[kMaxParentChildLinks];
  double fSeparation3D[kMaxParentChildLinks];
  double fEnergyRatio[kMaxParentChildLinks];
  double fPIDLinkType[kMaxParentChildLinks];
  double fOpeningAngle[kMaxParentChildLinks];
  double fTrackShowerLinkType[kMaxParentChildLinks];
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
  //Modes
  ////////////////////////////////////////
  bool fVisualisationMode;
  ////////////////////////////////////////
  //Algs
  ////////////////////////////////////////
  PandizzleAlg fPandizzleAlg;
  PandrizzleAlg fPandrizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  ctp::CTPHelper fConvTrackPID;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
  ivysaurus::IvysaurusEvaluator fIvysaurusEvaluator;
  TrueEnergyCalc fTrueEnergyCalc;
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
  fVisualisationMode(pset.get<bool>("VisualisationMode")),
  fPandizzleAlg(pset.get<fhicl::ParameterSet>("PandizzleConfig")),
  fPandrizzleAlg(pset.get<fhicl::ParameterSet>("PandrizzleConfig")),
  fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
  fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"), fTrackModuleLabel, fShowerModuleLabel,
      fHitsModuleLabel, fWireModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fRecoModuleLabel),
  fIvysaurusEvaluator(pset.get<fhicl::ParameterSet>("IvysaurusEvaluator")),
  fTrueEnergyCalc(pset.get<fhicl::ParameterSet>("TrueEnergyCalc"))
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
    std::cout << "AAAAAA" << std::endl;
    Reset();
    std::cout << "BBB" << std::endl;
    fRun = evt.run();
    fSubRun = evt.subRun();
    fEvent = evt.event();
    fIsMC = !evt.isRealData();
    std::cout << "CCC" << std::endl;
    FillPandoraMaps(evt);
    std::cout << "DDD" << std::endl;
    GetEventInfo(evt);
    std::cout << "EEE" << std::endl;
    if (fIsMC) 
        GetTruthInfo(evt);
    std::cout << "FFF" << std::endl;
    FillVertexInfo(evt);
    std::cout << "GGG" << std::endl;
    FillPFParticleInfo(evt);
    std::cout << "HHH" << std::endl;
    FillHierarchyInfo(evt);
    std::cout << "III" << std::endl;
    FillParentChildLinkInfo(evt);
    std::cout << "JJJ" << std::endl;
    RunTrackSelection(evt);
    std::cout << "KKK" << std::endl;
    RunShowerSelection(evt);
    std::cout << "LLL" << std::endl;

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
    fPOTTree->Branch("POT", &fPOT);
    fPOTTree->Branch("Run", &fRun);
    fPOTTree->Branch("SubRun", &fSubRun);

    ///////////////////////////
    // CCNuSelection tree
    ///////////////////////////
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");

    ////////////////////////////
    // Event info
    ////////////////////////////
    fTree->Branch("Run", &fRun);
    fTree->Branch("SubRun", &fSubRun);
    fTree->Branch("Event", &fEvent);

    ////////////////////////////
    // True info
    ////////////////////////////
    fTree->Branch("NuPdg", &fNuPdg);
    fTree->Branch("BeamPdg", &fBeamPdg);
    fTree->Branch("NuTrackID", &fNuTrackID);
    fTree->Branch("NC", &fNC);
    fTree->Branch("Mode", &fMode);
    fTree->Branch("TargetZ", &fTargetZ);
    fTree->Branch("Enu", &fENu);
    fTree->Branch("TrueNuEnergyEDep", &fTrueNuEnergyEDep);
    fTree->Branch("NuX", &fNuX);
    fTree->Branch("NuY", &fNuY);
    fTree->Branch("NuZ", &fNuZ);
    fTree->Branch("IsMC", &fIsMC);
    fTree->Branch("T0", &fT0);
    fTree->Branch("Q2", &fQ2);
    fTree->Branch("W", &fW);
    fTree->Branch("X", &fX);
    fTree->Branch("Y", &fY);
    fTree->Branch("NuMomX", &fNuMomX);
    fTree->Branch("NuMomY", &fNuMomY);
    fTree->Branch("NuMomZ", &fNuMomZ);
    fTree->Branch("NuMomT", &fNuMomT);
    fTree->Branch("NuT", &fNuT);
    fTree->Branch("LepPDG", &fLepPDG);
    fTree->Branch("LepEnergy", &fLepEnergy);
    fTree->Branch("TrueLepEnergyEDep", &fTrueLepEnergyEDep);
    fTree->Branch("MomLepX", &fMomLepX);
    fTree->Branch("MomLepY", &fMomLepY);
    fTree->Branch("MomLepZ", &fMomLepZ);
    fTree->Branch("MomLepT", &fMomLepT);
    fTree->Branch("LepEndX", &fLepEndX);
    fTree->Branch("LepEndY", &fLepEndY);
    fTree->Branch("LepEndZ", &fLepEndZ);
    fTree->Branch("LepEndT", &fLepEndT);
    fTree->Branch("LepNuAngle", &fLepNuAngle);
    fTree->Branch("LepHasMichelDecay", &fLepHasMichelDecay);
    fTree->Branch("NuMomTranMag", &fNuMomTranMag);
    fTree->Branch("TargNuclMomTranMag", &fTargNuclMomTranMag);
    fTree->Branch("InitalMomTranMag", &fInitalMomTranMag);
    fTree->Branch("LepMomTranMag", &fLepMomTranMag);
    fTree->Branch("NuclRemMomTranMag", &fNuclRemMomTranMag);
    fTree->Branch("FinalMomTranMagNoLepNoRem", &fFinalMomTranMagNoLepNoRem);
    fTree->Branch("FinalMomTranMagNoLepWithRem", &fFinalMomTranMagNoLepWithRem);
    fTree->Branch("FinalMomTranMagWithLepNoRem", &fFinalMomTranMagWithLepNoRem);
    fTree->Branch("FinalMomTranMagWithLepWithRem", &fFinalMomTranMagWithLepWithRem);
    fTree->Branch("TrueEnergyDepX_Inner", &fTrueEnergyDepX_Inner);
    fTree->Branch("TrueEnergyDepY_Inner", &fTrueEnergyDepY_Inner);
    fTree->Branch("TrueEnergyDepZ_Inner", &fTrueEnergyDepZ_Inner);
    fTree->Branch("IsEdepInAPA", &fIsEdepInAPA);
    fTree->Branch("APALowY", &fAPALowY);
    fTree->Branch("APAHighY", &fAPAHighY);
    fTree->Branch("APALowZ", &fAPALowZ);
    fTree->Branch("APAHighZ", &fAPAHighZ);

    ////////////////////////////
    // Event-level reco info 
    ////////////////////////////
    fTree->Branch("NRecoPFPs", &fNRecoPFPs);
    fTree->Branch("RecoNuVtxX", &fRecoNuVtxX);
    fTree->Branch("RecoNuVtxY", &fRecoNuVtxY);
    fTree->Branch("RecoNuVtxZ", &fRecoNuVtxZ);
    fTree->Branch("RecoNuVtxNShowers", &fRecoNuVtxNShowers);
    fTree->Branch("RecoNuVtxNTracks", &fRecoNuVtxNTracks);
    fTree->Branch("RecoNuVtxNChildren", &fRecoNuVtxNChildren);
    fTree->Branch("RecoEventCharge", &fRecoEventCharge);
    fTree->Branch("NumuRecoENu", fNumuRecoENu, "NumuRecoENu[NRecoPFPs]/D");
    fTree->Branch("NumuRecoMomLep", fNumuRecoMomLep, "NumuRecoMomLep[NRecoPFPs]/D");
    fTree->Branch("NumuRecoEHad", fNumuRecoEHad, "fNumuRecoEHad[NRecoPFPs]/D");
    fTree->Branch("NueRecoENu", fNueRecoENu, "fNueRecoENu[NRecoPFPs]/D");
    fTree->Branch("NueRecoMomLep", fNueRecoMomLep, "fNueRecoMomLep[NRecoPFPs]/D");
    fTree->Branch("NueRecoEHad", fNueRecoEHad, "fNueRecoEHad[NRecoPFPs]/D");

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
    // True
    fTree->Branch("RecoPFPTruePDG", fRecoPFPTruePDG, "RecoPFPTruePDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueTrackID", fRecoPFPTrueTrackID, "RecoPFPTrueTrackID[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTruePrimary", fRecoPFPTruePrimary,"RecoPFPTruePrimary[NRecoPFPs]/O");
    fTree->Branch("RecoPFPTrueGeneration", fRecoPFPTrueGeneration, "RecoPFPTrueGeneration[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueParentTrackID", fRecoPFPTrueParentTrackID, "RecoPFPTrueParentTrackID[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueParentPDG", fRecoPFPTrueParentPDG, "RecoPFPTrueParentPDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueVisibleGeneration", fRecoPFPTrueVisibleGeneration, "RecoPFPTrueVisibleGeneration[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueVisibleParentTrackID", fRecoPFPTrueVisibleParentTrackID, "RecoPFPTrueVisibleParentTrackID[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueVisibleParentPDG", fRecoPFPTrueVisibleParentPDG, "RecoPFPTrueVisibleParentPDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrueEnergyEDep", fRecoPFPTrueEnergyEDep, "RecoPFPTrueEnergyEDep[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomX", fRecoPFPTrueMomX,"RecoPFPTrueMomX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomY", fRecoPFPTrueMomY,"RecoPFPTrueMomY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomZ", fRecoPFPTrueMomZ,"RecoPFPTrueMomZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueMomT", fRecoPFPTrueMomT,"RecoPFPTrueMomT[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartX", fRecoPFPTrueStartX,"RecoPFPTrueStartX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartY", fRecoPFPTrueStartY,"RecoPFPTrueStartY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueStartZ", fRecoPFPTrueStartZ,"RecoPFPTrueStartZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndX", fRecoPFPTrueEndX,"RecoPFPTrueEndX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndY", fRecoPFPTrueEndY,"RecoPFPTrueEndY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrueEndZ", fRecoPFPTrueEndZ,"RecoPFPTrueEndZ[NRecoPFPs]/D");
    // Reco
    fTree->Branch("RecoPFPSelf", fRecoPFPSelf, "RecoPFPSelf[NRecoPFPs]/I");
    fTree->Branch("RecoPFPIsPrimary", fRecoPFPIsPrimary, "RecoPFPIsPrimary[NRecoPFPs]/O");
    fTree->Branch("RecoPFPRecoGeneration", fRecoPFPRecoGeneration, "RecoPFPRecoGeneration[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoParentSelf", fRecoPFPRecoParentSelf, "RecoPFPRecoParentSelf[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoParentPDG", fRecoPFPRecoParentPDG, "RecoPFPRecoParentPDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrackShowerScore", fRecoPFPTrackShowerScore, "RecoPFPTrackShowerScore[NRecoPFPs]/D");
    fTree->Branch("RecoPFPTrackShowerPDG", fRecoPFPTrackShowerPDG, "RecoPFPTrackShowerPDG[NRecoPFPs]/I");
    fTree->Branch("RecoPFPSpacepointX", &fRecoPFPSpacepointX);
    fTree->Branch("RecoPFPSpacepointY", &fRecoPFPSpacepointY);
    fTree->Branch("RecoPFPSpacepointZ", &fRecoPFPSpacepointZ);
    fTree->Branch("RecoPFPRecoNHits", fRecoPFPRecoNHits,"RecoPFPRecoNHits[NRecoPFPs]/I");
    fTree->Branch("RecoPFPShowerFitSuccess", fRecoPFPShowerFitSuccess, "RecoPFPShowerFitSuccess[NRecoPFPs]/I");
    fTree->Branch("RecoPFPTrackFitSuccess", fRecoPFPTrackFitSuccess, "RecoPFPTrackFitSuccess[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoCompleteness", fRecoPFPRecoCompleteness,"RecoPFPRecoCompleteness[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoHitPurity", fRecoPFPRecoHitPurity,"RecoPFPRecoHitPurity[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoCharge", fRecoPFPRecoCharge,"RecoPFPRecoCharge[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexX", fRecoPFPRecoVertexX,"RecoPFPRecoVertexX[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexY", fRecoPFPRecoVertexY,"RecoPFPRecoVertexY[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoVertexZ", fRecoPFPRecoVertexZ,"RecoPFPRecoVertexZ[NRecoPFPs]/D");
    fTree->Branch("RecoPFPRecoNChildPFP", fRecoPFPRecoNChildPFP,"RecoPFPRecoNChildPFP[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoNChildTrackPFP", fRecoPFPRecoNChildTrackPFP,"RecoPFPRecoNChildTrackPFP[NRecoPFPs]/I");
    fTree->Branch("RecoPFPRecoNChildShowerPFP", fRecoPFPRecoNChildShowerPFP,"RecoPFPRecoNChildShowerPFP[NRecoPFPs]/I");
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
    fTree->Branch("RecoPFPIvysaurusParticleType", fRecoPFPIvysaurusParticleType, "RecoPFPIvysaurusParticleType[NRecoPFPs]/I");

    ////////////////////////////
    // Track Info
    ////////////////////////////
    // Reco
    fTree->Branch("RecoTrackRecoStartX", fRecoTrackRecoStartX, "RecoTrackRecoStartX[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoStartY", fRecoTrackRecoStartY, "RecoTrackRecoStartY[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoStartZ", fRecoTrackRecoStartZ, "RecoTrackRecoStartZ[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndX", fRecoTrackRecoEndX, "RecoTrackRecoEndX[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndY", fRecoTrackRecoEndY, "RecoTrackRecoEndY[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndZ", fRecoTrackRecoEndZ, "RecoTrackRecoEndZ[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoUpstreamX", fRecoTrackRecoUpstreamX, "RecoTrackRecoUpstreamX[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoUpstreamY", fRecoTrackRecoUpstreamY, "RecoTrackRecoUpstreamY[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoUpstreamZ", fRecoTrackRecoUpstreamZ, "RecoTrackRecoUpstreamZ[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoDownstreamX", fRecoTrackRecoDownstreamX, "RecoTrackRecoDownstreamX[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoDownstreamY", fRecoTrackRecoDownstreamY, "RecoTrackRecoDownstreamY[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoDownstreamZ", fRecoTrackRecoDownstreamZ, "RecoTrackRecoDownstreamZ[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexX", fRecoTrackRecoEndClosestToVertexX, "RecoTrackRecoEndClosestToVertexX[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexY", fRecoTrackRecoEndClosestToVertexY, "RecoTrackRecoEndClosestToVertexY[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexZ", fRecoTrackRecoEndClosestToVertexZ, "RecoTrackRecoEndClosestToVertexZ[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoLength", fRecoTrackRecoLength, "RecoTrackRecoLength[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoContained", fRecoTrackRecoContained, "RecoTrackRecoContained[NRecoPFPs]/I");
    fTree->Branch("RecoTrackRecoMomMethod", fRecoTrackRecoMomMethod, "RecoTrackRecoMomMethod[NRecoPFPs]/I");
    fTree->Branch("RecoTrackRecoMomMCS", fRecoTrackRecoMomMCS, "RecoTrackRecoMomMCS[NRecoPFPs]/D");
    fTree->Branch("RecoTrackRecoMomRange", fRecoTrackRecoMomRange, "RecoTrackRecoMomRange[NRecoPFPs]/D");
    // Pandizzle
    fTree->Branch("RecoTrackMichelNHits", fRecoTrackMichelNHits, "RecoTrackMichelNHits[NRecoPFPs]/D");
    fTree->Branch("RecoTrackMichelElectronMVA", fRecoTrackMichelElectronMVA, "RecoTrackMichelElectronMVA[NRecoPFPs]/D");
    fTree->Branch("RecoTrackMichelRecoEnergyPlane2", fRecoTrackMichelRecoEnergyPlane2, "RecoTrackMichelRecoEnergyPlane2[NRecoPFPs]/D");
    fTree->Branch("RecoTrackDeflecAngleSD", fRecoTrackDeflecAngleSD, "RecoTrackDeflecAngleSD[NRecoPFPs]/D");
    fTree->Branch("RecoTrackLength", fRecoTrackLength, "RecoTrackLength[NRecoPFPs]/D");
    fTree->Branch("RecoTrackEvalRatio", fRecoTrackEvalRatio, "RecoTrackEvalRatio[NRecoPFPs]/D");
    fTree->Branch("RecoTrackConcentration", fRecoTrackConcentration, "RecoTrackConcentration[NRecoPFPs]/D");
    fTree->Branch("RecoTrackCoreHaloRatio", fRecoTrackCoreHaloRatio, "RecoTrackCoreHaloRatio[NRecoPFPs]/D");
    fTree->Branch("RecoTrackConicalness", fRecoTrackConicalness, "RecoTrackConicalness[NRecoPFPs]/D");
    fTree->Branch("RecoTrackdEdxStart", fRecoTrackdEdxStart, "RecoTrackdEdxStart[NRecoPFPs]/D");
    fTree->Branch("RecoTrackdEdxEnd", fRecoTrackdEdxEnd, "RecoTrackdEdxEnd[NRecoPFPs]/D");
    fTree->Branch("RecoTrackdEdxEndRatio", fRecoTrackdEdxEndRatio, "RecoTrackdEdxEndRatio[NRecoPFPs]/D");
    fTree->Branch("RecoTrackPandizzleVar", fRecoTrackPandizzleVar, "RecoTrackPandizzleVar[NRecoPFPs]/D");
    // Selected info
    fTree->Branch("SelTrackPandizzleSelf", &fSelTrackPandizzleSelf);
    fTree->Branch("SelTrackPandizzleIndex", &fSelTrackPandizzleIndex);
    fTree->Branch("SelTrackDeepPanSelf", &fSelTrackDeepPanSelf);
    fTree->Branch("SelTrackDeepPanIndex", &fSelTrackDeepPanIndex);
    fTree->Branch("SelTrackIvysaurusSelf", &fSelTrackIvysaurusSelf);
    fTree->Branch("SelTrackIvysaurusIndex", &fSelTrackIvysaurusIndex);
    fTree->Branch("SelTrackLongestLengthSelf", &fSelTrackLongestLengthSelf);
    fTree->Branch("SelTrackLongestLengthIndex", &fSelTrackLongestLengthIndex);
    ///////////////////////////
    // Shower Info
    ///////////////////////////
    // Reco
    fTree->Branch("RecoShowerRecoStartX", fRecoShowerRecoStartX, "RecoShowerRecoStartX[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoStartY", fRecoShowerRecoStartY, "RecoShowerRecoStartY[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoStartZ", fRecoShowerRecoStartZ, "RecoShowerRecoStartZ[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoDirX", fRecoShowerRecoDirX, "RecoShowerRecoDirX[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoDirY", fRecoShowerRecoDirY, "RecoShowerRecoDirY[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoDirZ", fRecoShowerRecoDirZ, "RecoShowerRecoDirZ[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoLength", &fRecoShowerRecoLength, "RecoShowerRecoLength[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoOpeningAngle", &fRecoShowerRecoOpeningAngle, "RecoShowerRecoOpeningAngle[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecodEdx", fRecoShowerRecodEdx, "RecoShowerRecodEdx[NRecoPFPs][3]/D");
    fTree->Branch("RecoShowerRecoBestPlane", &fRecoShowerRecoBestPlane, "RecoShowerRecoBestPlane[NRecoPFPs]/I");
    fTree->Branch("RecoShowerRecoMom", fRecoShowerRecoMom, "RecoShowerRecoMom[NRecoPFPs]/D");
    fTree->Branch("RecoShowerRecoEnergy", fRecoShowerRecoEnergy, "RecoShowerRecoEnergy[NRecoPFPs][3]/D");
    // Pandrizzle
    fTree->Branch("RecoShowerPandrizzleConnectionBDTScore", &fRecoShowerPandrizzleConnectionBDTScore, "RecoShowerPandrizzleConnectionBDTScore[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzlePathwayLengthMin", &fRecoShowerPandrizzlePathwayLengthMin, "RecoShowerPandrizzlePathwayLengthMin[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D", &fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D, 
                  "RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxNPostShowerStartHits", &fRecoShowerPandrizzleMaxNPostShowerStartHits, "RecoShowerPandrizzleMaxNPostShowerStartHits[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartScatterAngle", &fRecoShowerPandrizzleMaxPostShowerStartScatterAngle, "RecoShowerPandrizzleMaxPostShowerStartScatterAngle[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry, 
                  "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry, 
                  "RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance, 
                  "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius", &fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius, 
                  "RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartOpeningAngle", &fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle, "RecoShowerPandrizzleMaxPostShowerStartOpeningAngle[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxFoundHitRatio", &fRecoShowerPandrizzleMaxFoundHitRatio, "RecoShowerPandrizzleMaxFoundHitRatio[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMaxInitialGapSize", &fRecoShowerPandrizzleMaxInitialGapSize, "RecoShowerPandrizzleMaxInitialGapSize[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleMinLargestProjectedGapSize", &fRecoShowerPandrizzleMinLargestProjectedGapSize, "RecoShowerPandrizzleMinLargestProjectedGapSize[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleNViewsWithAmbiguousHits", &fRecoShowerPandrizzleNViewsWithAmbiguousHits, "RecoShowerPandrizzleNViewsWithAmbiguousHits[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy", &fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy, "RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleEvalRatio", &fRecoShowerPandrizzleEvalRatio, "RecoShowerPandrizzleEvalRatio[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleConcentration", &fRecoShowerPandrizzleConcentration, "RecoShowerPandrizzleConcentration[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleCoreHaloRatio", &fRecoShowerPandrizzleCoreHaloRatio, "RecoShowerPandrizzleCoreHaloRatio[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleConicalness", &fRecoShowerPandrizzleConicalness, "RecoShowerPandrizzleConicalness[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzledEdxBestPlane", &fRecoShowerPandrizzledEdxBestPlane, "RecoShowerPandrizzledEdxBestPlane[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleDisplacement", &fRecoShowerPandrizzleDisplacement, "RecoShowerPandrizzleDisplacement[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleDCA", &fRecoShowerPandrizzleDCA, "RecoShowerPandrizzleDCA[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleWideness", &fRecoShowerPandrizzleWideness, "RecoShowerPandrizzleWideness[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleEnergyDensity", &fRecoShowerPandrizzleEnergyDensity, "RecoShowerPandrizzleEnergyDensity[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleBDTMethod", &fRecoShowerPandrizzleBDTMethod, "RecoShowerPandrizzleBDTMethod[NRecoPFPs]/D");
    fTree->Branch("RecoShowerEnhancedPandrizzleScore", &fRecoShowerEnhancedPandrizzleScore, "RecoShowerEnhancedPandrizzleScore[NRecoPFPs]/D"); 
    fTree->Branch("RecoShowerBackupPandrizzleScore", &fRecoShowerBackupPandrizzleScore, "RecoShowerBackupPandrizzleScore[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleIsFilled", &fRecoShowerPandrizzleIsFilled, "RecoShowerPandrizzleIsFilled[NRecoPFPs]/O");
    fTree->Branch("RecoShowerPandrizzleModularPathwayLength", &fRecoShowerPandrizzleModularPathwayLength, "RecoShowerPandrizzleModularPathwayLength[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance", &fRecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance, "RecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance[NRecoPFPs]/D");
    fTree->Branch("RecoShowerPandrizzleModularMaxNShowerHits", &fRecoShowerPandrizzleModularMaxNShowerHits, "RecoShowerPandrizzleModularMaxNShowerHits[NRecoPFPs]/D");

    // Selected info
    fTree->Branch("SelShowerPandrizzleSelf", &fSelShowerPandrizzleSelf);
    fTree->Branch("SelShowerPandrizzleIndex", &fSelShowerPandrizzleIndex);
    fTree->Branch("SelShowerIvysaurusSelf", &fSelShowerIvysaurusSelf);
    fTree->Branch("SelShowerIvysaurusIndex", &fSelShowerIvysaurusIndex);
    fTree->Branch("SelShowerHighestEnergySelf", &fSelShowerHighestEnergySelf);
    fTree->Branch("SelShowerHighestEnergyIndex", &fSelShowerHighestEnergyIndex);

    ////////////////////////////////////////
    // Hierarchy Info
    ////////////////////////////////////////
    fTree->Branch("NParentChildLinks", &fNParentChildLinks);
    // Truth
    fTree->Branch("TrueParentChildLink", &fTrueParentChildLink, "TrueParentChildLink[NParentChildLinks]/O");
    // Parent information
    fTree->Branch("ParentTrackScore", &fParentTrackScore, "ParentTrackScore[NParentChildLinks]/D");
    fTree->Branch("ParentBraggVariable", &fParentBraggVariable, "ParentBraggVariable[NParentChildLinks]/D");
    fTree->Branch("ParentEndRegionNHits", &fParentEndRegionNHits, "ParentEndRegionNHits[NParentChildLinks]/D");
    fTree->Branch("ParentEndRegionNParticles", &fParentEndRegionNParticles, "ParentEndRegionNParticles[NParentChildLinks]/D");
    fTree->Branch("ParentEndRegionRToWall", &fParentEndRegionRToWall, "ParentEndRegionRToWall[NParentChildLinks]/D");
    // Edge information
    fTree->Branch("ParentPFPIndex", &fParentPFPIndex, "ParentPFPIndex[NParentChildLinks]/D");
    fTree->Branch("ChildPFPIndex", &fChildPFPIndex, "ChildPFPIndex[NParentChildLinks]/D");
    fTree->Branch("VertexSeparation", &fVertexSeparation, "VertexSeparation[NParentChildLinks]/D");
    fTree->Branch("SeparationU", &fSeparationU, "SeparationU[NParentChildLinks]/D");
    fTree->Branch("SeparationV", &fSeparationV, "SeparationV[NParentChildLinks]/D");
    fTree->Branch("SeparationW", &fSeparationW, "SeparationW[NParentChildLinks]/D");
    fTree->Branch("Separation3D", &fSeparation3D, "Separation3D[NParentChildLinks]/D");
    fTree->Branch("EnergyRatio", &fEnergyRatio, "EnergyRatio[NParentChildLinks]/D");
    fTree->Branch("PIDLinkType", &fPIDLinkType, "PIDLinkType[NParentChildLinks]/D");
    fTree->Branch("OpeningAngle", &fOpeningAngle, "OpeningAngle[NParentChildLinks]/D");
    fTree->Branch("TrackShowerLinkType", &fTrackShowerLinkType, "TrackShowerLinkType[NParentChildLinks]/D");
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

  if (sr.getByLabel(fPOTModuleLabel, potListHandle))
  {
    fPOT = potListHandle->totpot;
  }
  else
  {
    fPOT = 0.;
  }

  if (fPOTTree) fPOTTree->Fill();
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::endJob()
{
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::Reset()
{
    fTrueEnergyCalc.Reset();

    ////////////////////////////////////////
    // Pandora maps
    ////////////////////////////////////////
    fMCParticleMap.clear();
    fPFPMap.clear();
    fReconstructedMCParticles.clear();
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
    fNuTrackID = kDefInt;
    fNC = kDefInt;    
    fMode = kDefInt; 
    fTargetZ = kDefInt;
    fQ2 = kDefDoub; 
    fENu = kDefDoub; 
    fTrueNuEnergyEDep = kDefDoub;
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
    fLepPDG = kDefInt;
    fLepEnergy = kDefDoub;
    fTrueLepEnergyEDep = kDefDoub;
    fMomLepX = kDefDoub;
    fMomLepY = kDefDoub;
    fMomLepZ = kDefDoub;
    fMomLepT = kDefDoub;
    fLepEndX = kDefDoub;
    fLepEndY = kDefDoub;
    fLepEndZ = kDefDoub;
    fLepEndT = kDefDoub;
    fLepNuAngle = kDefDoub;
    fLepHasMichelDecay = 0;
    fNuMomTranMag = kDefDoub;
    fTargNuclMomTranMag = kDefDoub;
    fInitalMomTranMag = kDefDoub;
    fLepMomTranMag = kDefDoub;
    fNuclRemMomTranMag = kDefDoub;
    fFinalMomTranMagNoLepNoRem = kDefDoub;
    fFinalMomTranMagNoLepWithRem = kDefDoub;
    fFinalMomTranMagWithLepNoRem = kDefDoub;
    fFinalMomTranMagWithLepWithRem = kDefDoub;
    fTrueEnergyDepX_Inner.clear();
    fTrueEnergyDepY_Inner.clear();
    fTrueEnergyDepZ_Inner.clear();
    fIsEdepInAPA.clear();
    fAPALowY.clear();
    fAPAHighY.clear();
    fAPALowZ.clear();
    fAPAHighZ.clear();

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
        // Event-level reco info 
        ////////////////////////////
        fNumuRecoENu[i] = kDefDoub;
        fNumuRecoMomLep[i] = kDefDoub;
        fNumuRecoEHad[i] = kDefDoub;
        fNueRecoENu[i] = kDefDoub;
        fNueRecoMomLep[i] = kDefDoub;
        fNueRecoEHad[i] = kDefDoub;

        ////////////////////////////  
        // PFParticle stuff
        ////////////////////////////  
        // Truth
        fRecoPFPTruePDG[i] = kDefInt;
        fRecoPFPTrueTrackID[i] = kDefInt;
        fRecoPFPTruePrimary[i] = false;
        fRecoPFPTrueGeneration[i] = kDefInt;
        fRecoPFPTrueParentTrackID[i] = kDefInt;
        fRecoPFPTrueParentPDG[i] = kDefInt;
        fRecoPFPTrueVisibleGeneration[i] = kDefInt;
        fRecoPFPTrueVisibleParentTrackID[i] = kDefInt;
        fRecoPFPTrueVisibleParentPDG[i] = kDefInt;
        fRecoPFPTrueEnergyEDep[i] = kDefDoub;
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
        fRecoPFPRecoGeneration[i] = kDefInt;
        fRecoPFPRecoParentSelf[i] = kDefInt;
        fRecoPFPRecoParentPDG[i] = kDefInt;
        fRecoPFPTrackShowerScore[i] = kDefDoub;
        fRecoPFPTrackShowerPDG[i] = kDefInt;
        fRecoPFPRecoNHits[i] = kDefInt;
        fRecoPFPSpacepointX.clear();
        fRecoPFPSpacepointY.clear();
        fRecoPFPSpacepointZ.clear();
        fRecoPFPShowerFitSuccess[i] = 0;
        fRecoPFPTrackFitSuccess[i] = 0;
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
        fRecoTrackRecoMomRange[i] = kDefDoub;
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
        fRecoShowerPandrizzleModularPathwayLength[i] = kDefDoub;
        fRecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance[i] = kDefDoub;
        fRecoShowerPandrizzleModularMaxNShowerHits[i] = kDefDoub;
        fRecoShowerPandrizzleBDTMethod[i] = kDefDoub;
        fRecoShowerEnhancedPandrizzleScore[i] = kDefDoub;
        fRecoShowerBackupPandrizzleScore[i] = kDefDoub;
        fRecoShowerPandrizzleIsFilled[i] = 0;
    }
    // Selection stuff
    fSelTrackPandizzleSelf = -1;
    fSelTrackPandizzleIndex = -1;
    fSelTrackDeepPanSelf = -1;
    fSelTrackDeepPanIndex = -1;
    fSelTrackIvysaurusSelf = -1;
    fSelTrackIvysaurusIndex = -1;
    fSelTrackLongestLengthSelf = -1;
    fSelTrackLongestLengthIndex = -1;
    fSelShowerPandrizzleSelf = -1;
    fSelShowerPandrizzleIndex = -1;
    fSelShowerIvysaurusSelf = -1;
    fSelShowerIvysaurusIndex = -1;
    fSelShowerHighestEnergySelf = -1;
    fSelShowerHighestEnergyIndex = -1;


    ////////////////////////////////////////
    // Hierarchy Info
    ////////////////////////////////////////
    // Reco
    fNParentChildLinks = 0;

    for (int i = 0; i < kMaxParentChildLinks; i++)
    {
        // Truth
        fTrueParentChildLink[i] = false;
        // Parent information
        fParentTrackScore[i] = kDefDoub;
        fParentBraggVariable[i] = kDefDoub;
        fParentEndRegionNHits[i] = kDefDoub;
        fParentEndRegionNParticles[i] = kDefDoub;
        fParentEndRegionRToWall[i] = kDefDoub;
        // Edge information
        fParentPFPIndex[i] = kDefDoub;
        fChildPFPIndex[i] = kDefDoub;
        fVertexSeparation[i] = kDefDoub;
        fSeparationU[i] = kDefDoub;
        fSeparationV[i] = kDefDoub;
        fSeparationW[i] = kDefDoub;
        fSeparation3D[i] = kDefDoub;
        fEnergyRatio[i] = kDefDoub;
        fPIDLinkType[i] = kDefDoub;
        fOpeningAngle[i] = kDefDoub;
        fTrackShowerLinkType[i] = kDefDoub;
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
    fTrueEnergyCalc.InitialiseCalc(evt);

    /*
    // Investigating stuff...
    //////////////////////////////////////////////////////////////////////////////////
    // Get the 'inner' sim energy deposit positions
    art::Handle<std::vector<sim::SimEnergyDeposit>> simEDepHandle_inner;
    std::vector<art::Ptr<sim::SimEnergyDeposit>> simEDepVector_inner;

    if (evt.getByLabel("largeant:LArG4DetectorServicevolTPCActiveInner", simEDepHandle_inner))
        art::fill_ptr_vector(simEDepVector_inner, simEDepHandle_inner);

    lar_pandora::LArDriftVolumeList driftVolumeList; // list of 'drift volumes' (i think there's a bug.. there's a lot of repitition in this map)
    lar_pandora::LArDriftVolumeMap tpcVolumeMap; // map of 'tpc volumes' to 'tpc volume list!' i.e. apa volumes
    bool useActiveBoundingBox = true; //true;

    lar_pandora::LArPandoraGeometry::LoadGeometry(driftVolumeList, tpcVolumeMap, useActiveBoundingBox);

    std::vector<float> yGapBoundaries, zGapBoundaries;

    for (const auto &driftVol : driftVolumeList)
    {
        for (const auto &apaVol : driftVol.GetTpcVolumeList())
        {
            const float yLow = apaVol.GetCenterY() - (apaVol.GetWidthY() * 0.5);
            const float yHigh = apaVol.GetCenterY() + (apaVol.GetWidthY() * 0.5);
            const float zLow = apaVol.GetCenterZ() - (apaVol.GetWidthZ() * 0.5);
            const float zHigh = apaVol.GetCenterZ() + (apaVol.GetWidthZ() * 0.5);

            fAPALowY.push_back(yLow);
            fAPAHighY.push_back(yHigh);
            fAPALowZ.push_back(zLow);
            fAPAHighZ.push_back(zHigh);
        }
    }

    for (const art::Ptr<sim::SimEnergyDeposit> &simEDep : simEDepVector_inner)
    {
        const float x = simEDep->X();
        const float y = simEDep->Y();
        const float z = simEDep->Z();

        fTrueEnergyDepX_Inner.push_back(x);
        fTrueEnergyDepY_Inner.push_back(y);
        fTrueEnergyDepZ_Inner.push_back(z);

        int isInAPAVol = 0;

        for (unsigned int iAPA = 0; iAPA < fAPALowY.size(); iAPA++)
        {
            const float apaLowY = fAPALowY.at(iAPA);
            const float apaHighY = fAPAHighY.at(iAPA);
            const float apaLowZ = fAPALowZ.at(iAPA);
            const float apaHighZ = fAPAHighZ.at(iAPA);

            if ((y > apaLowY) && (y < apaHighY) && (z > apaLowZ) && (z < apaHighZ))
            {
                isInAPAVol = 1;
                break;
            }
        }

        fIsEdepInAPA.push_back(isInAPAVol);
    }
    */
    //////////////////////////////////////////////////////////////////////////////////


    const std::vector<art::Ptr<simb::MCTruth>> mcTruths = dune_ana::DUNEAnaEventUtils::GetMCTruths(evt, fNuGenModuleLabel);

    if (mcTruths.empty())
        return;

    art::Ptr<simb::MCTruth> mcTruth = mcTruths.at(0);

    if (mcTruth->Origin() != simb::kBeamNeutrino)
        return;

    // Get the original neutrino flavour (before osc.)
    art::Handle<std::vector<simb::MCFlux>> mcFluxListHandle;
    std::vector<art::Ptr<simb::MCFlux>> mcFlux;

    if (evt.getByLabel(fNuGenModuleLabel, mcFluxListHandle))
        art::fill_ptr_vector(mcFlux, mcFluxListHandle);

    if (mcFluxListHandle.isValid()) 
        fBeamPdg  = mcFlux[0]->fntype;

    // Neutrino
    const simb::MCNeutrino &mcNeutrino = mcTruth->GetNeutrino();
    fNuTrackID = mcNeutrino.Nu().TrackId();
    fNuPdg = mcNeutrino.Nu().PdgCode();
    fNC = mcNeutrino.CCNC();
    fMode = mcNeutrino.Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    fTargetZ = mcNeutrino.Target()%100000000/10000;
    fENu = mcNeutrino.Nu().E();
    fTrueNuEnergyEDep = fTrueEnergyCalc.GetTrueNuEnergy();
    fQ2 = mcNeutrino.QSqr();
    fW = mcNeutrino.W();
    fX = mcNeutrino.X();
    fY = mcNeutrino.Y();
    fNuX = mcNeutrino.Nu().Vx();
    fNuY = mcNeutrino.Nu().Vy();
    fNuZ = mcNeutrino.Nu().Vz();
    fNuT = mcNeutrino.Nu().T();
    fNuMomX = mcNeutrino.Nu().Momentum().X();
    fNuMomY = mcNeutrino.Nu().Momentum().Y();
    fNuMomZ = mcNeutrino.Nu().Momentum().Z();
    fNuMomT = mcNeutrino.Nu().Momentum().T();

    // Leading lepton
    const simb::MCParticle &mcLepton = mcNeutrino.Lepton();
    fLepPDG = mcLepton.PdgCode();
    fLepEnergy = mcLepton.E();
    fTrueLepEnergyEDep = fTrueEnergyCalc.GetTrueParticleEnergy(mcLepton.TrackId());
    fMomLepX = mcLepton.Momentum().X();
    fMomLepY = mcLepton.Momentum().Y();
    fMomLepZ = mcLepton.Momentum().Z();
    fMomLepT = mcLepton.Momentum().T();
    fLepEndX = mcLepton.EndPosition().X();
    fLepEndY = mcLepton.EndPosition().Y();
    fLepEndZ = mcLepton.EndPosition().Z();
    fLepEndY = mcLepton.EndPosition().T();
    fLepNuAngle = mcNeutrino.Nu().Momentum().Vect().Angle(mcLepton.Momentum().Vect());

    // Does leading lepton Michel decay?
    bool hasNumu = false, hasNue = false, hasElectron = false;

    if (std::abs(fLepPDG) == 13)
    {
        for (int i = 0; i < mcLepton.NumberDaughters(); ++i)
        {
            const int childTrackId = mcLepton.Daughter(i);

            if (fMCParticleMap.find(childTrackId) == fMCParticleMap.end())
                continue;

            const art::Ptr<simb::MCParticle> &childMCParticle = fMCParticleMap.at(childTrackId);
            const int childPDG = std::abs(childMCParticle->PdgCode());

            if (childPDG == 11)
                hasElectron = true;
            else if (childPDG == 12)
                hasNue = true;
            else if (childPDG == 14)
                hasNumu = true;
        }
    }

    fLepHasMichelDecay = (hasNumu && hasNue && hasElectron) ? 1 : 0;

    ///////////////////////////////////////////////////////////////////////////////////////////
    // DOM - do some transverse momentum stuff
    TVector3 beam_axis(0, 0, 1);
    beam_axis.RotateX(-0.101);
    TVector3 nu_mom_vect = mcNeutrino.Nu().Momentum().Vect();
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

    int pfpIndex = -1;
    fNRecoPFPs = 0;

    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        if ((std::fabs(pfp->PdgCode()) == 12) || (std::fabs(pfp->PdgCode()) == 14) || (std::fabs(pfp->PdgCode()) == 16))
            continue;

        pfpIndex++;
        fNRecoPFPs++;

        if (pfpIndex == kMaxPFParticles)
            break;

        fRecoPFPSelf[pfpIndex] = pfp->Self();

        for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
        {
            if (nuChild->Self() == pfp->Self())
            {
                fRecoPFPIsPrimary[pfpIndex] = true;
                break;
            }
        }

        fRecoPFPTrackShowerPDG[pfpIndex] = pfp->PdgCode();

        try
        {
            const art::Ptr<larpandoraobj::PFParticleMetadata> metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfp, evt, fRecoModuleLabel);

            if (metadata->GetPropertiesMap().find("TrackScore") != metadata->GetPropertiesMap().end())
                fRecoPFPTrackShowerScore[pfpIndex] = metadata->GetPropertiesMap().at("TrackScore");
        }
        catch (...) {}

        const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fRecoModuleLabel);
        fRecoPFPRecoNHits[pfpIndex] = pfpHits.size();


        // If visualise then drop in spacepoint information
        if (fVisualisationMode)
        {
            const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, fRecoModuleLabel);

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
        }

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
        fRecoPFPRecoCharge[pfpIndex]  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, pfpHits);

        try
        {
            art::Ptr<recob::Vertex> recoVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, fRecoModuleLabel);

            fRecoPFPRecoVertexX[pfpIndex] = recoVertex->position().X();
            fRecoPFPRecoVertexY[pfpIndex] = recoVertex->position().Y();
            fRecoPFPRecoVertexZ[pfpIndex] = recoVertex->position().Z();
        }
        catch(...)
        {
        }

        // Child particle info
        FillChildPFPInformation(pfp, evt, fRecoPFPRecoNChildPFP[pfpIndex], fRecoPFPRecoNChildTrackPFP[pfpIndex], fRecoPFPRecoNChildShowerPFP[pfpIndex]);

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
                fReconstructedMCParticles.push_back(fRecoPFPTrueTrackID[pfpIndex]);

                if (matched_mcparticle->Mother() == 0) 
                    fRecoPFPTruePrimary[pfpIndex] = true;
                else 
                    fRecoPFPTruePrimary[pfpIndex] = false;

                fRecoPFPTrueEnergyEDep[pfpIndex] = fTrueEnergyCalc.GetTrueParticleEnergy(matched_mcparticle->TrackId());
                fRecoPFPTrueMomX[pfpIndex] = matched_mcparticle->Momentum().X();
                fRecoPFPTrueMomY[pfpIndex] = matched_mcparticle->Momentum().Y();
                fRecoPFPTrueMomZ[pfpIndex] = matched_mcparticle->Momentum().Z();
                fRecoPFPTrueStartX[pfpIndex] = matched_mcparticle->Position(0).X();
                fRecoPFPTrueStartY[pfpIndex] = matched_mcparticle->Position(0).Y();
                fRecoPFPTrueStartZ[pfpIndex] = matched_mcparticle->Position(0).Z();
                fRecoPFPTrueEndX[pfpIndex] = matched_mcparticle->EndPosition().X();
                fRecoPFPTrueEndY[pfpIndex] = matched_mcparticle->EndPosition().Y();
                fRecoPFPTrueEndZ[pfpIndex] = matched_mcparticle->EndPosition().Z();
            }
        }

        // DeepPan stuff
        //ctp::CTPResult deepPanPIDResult = fConvTrackPID.RunConvolutionalTrackPID(pfp, evt);
        //if (deepPanPIDResult.IsValid())
        //{
        //fRecoPFPDeepPanMuVar[pfpIndex] = deepPanPIDResult.GetMuonScore();
        //fRecoPFPDeepPanPiVar[pfpIndex] = deepPanPIDResult.GetPionScore();
        //fRecoPFPDeepPanProtonVar[pfpIndex] = deepPanPIDResult.GetProtonScore();
        //}

        // Ivysaurus stuff
        ivysaurus::IvysaurusEvaluator::IvysaurusScores ivysaurusScores = fIvysaurusEvaluator.IvysaurusUseEvaluate(evt, pfp);
        fRecoPFPIvysaurusMuon[pfpIndex] = ivysaurusScores.m_muonScore;
        fRecoPFPIvysaurusProton[pfpIndex] =  ivysaurusScores.m_protonScore;
        fRecoPFPIvysaurusPion[pfpIndex] =  ivysaurusScores.m_pionScore;
        fRecoPFPIvysaurusElectron[pfpIndex] =  ivysaurusScores.m_electronScore;
        fRecoPFPIvysaurusPhoton[pfpIndex] =  ivysaurusScores.m_photonScore;
        fRecoPFPIvysaurusOther[pfpIndex] =  ivysaurusScores.m_otherScore;
        fRecoPFPIvysaurusParticleType[pfpIndex] =  ivysaurusScores.m_particleType;

        /*
        std::cout << "inside CCNuSel: " << std::endl;
        std::cout << "ivysaurusScores.m_muonScore: " << ivysaurusScores.m_muonScore << std::endl;
        std::cout << "ivysaurusScores.m_protonScore: " << ivysaurusScores.m_protonScore << std::endl;
        std::cout << "ivysaurusScores.m_pionScore: " << ivysaurusScores.m_pionScore << std::endl;
        std::cout << "ivysaurusScores.m_electronScore: " << ivysaurusScores.m_electronScore << std::endl;
        std::cout << "ivysaurusScores.m_photonScore: " << ivysaurusScores.m_photonScore << std::endl;
        std::cout << "ivysaurusScores.m_otherScore: " << ivysaurusScores.m_otherScore << std::endl; 
        */

        // Fill the track & shower information
        FillRecoTrackInfo(evt, pfp, pfpIndex);
        FillRecoShowerInfo(evt, pfp, pfpIndex);
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

void FDSelection::CCNuSelection::FillHierarchyInfo(art::Event const & evt)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
        return;

    std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);

    int pfpIndex = -1;

    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        if ((std::fabs(pfp->PdgCode()) == 12) || (std::fabs(pfp->PdgCode()) == 14) || (std::fabs(pfp->PdgCode()) == 16))
            continue;

        pfpIndex++;

        if (pfpIndex == kMaxPFParticles)
            break;

        SetRecoGenerationInfo(evt, pfp, fRecoPFPRecoGeneration[pfpIndex], fRecoPFPRecoParentSelf[pfpIndex], fRecoPFPRecoParentPDG[pfpIndex]);
        SetTrueGenerationInfo(evt, pfpIndex, false, fRecoPFPTrueGeneration[pfpIndex], fRecoPFPTrueParentTrackID[pfpIndex], fRecoPFPTrueParentPDG[pfpIndex]);
        SetTrueGenerationInfo(evt, pfpIndex, true, fRecoPFPTrueVisibleGeneration[pfpIndex], fRecoPFPTrueVisibleParentTrackID[pfpIndex], fRecoPFPTrueVisibleParentPDG[pfpIndex]);
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::SetRecoGenerationInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp, 
    int &recoGeneration, int &recoParentSelf, int &recoParentPDG)
{
    recoGeneration = lar_pandora::LArPandoraHelper::GetGeneration(fPFPMap, pfp);
    const int parentID = pfp->Parent();

    if ((recoGeneration == 1) || (fPFPMap.find(parentID) == fPFPMap.end()))
    {
        recoParentSelf = kDefInt;
        recoParentPDG = kDefInt;
        return;
    }

    recoParentSelf = fPFPMap.at(parentID)->Self();
    recoParentPDG = kDefInt; // for now.. else I'll have to implement a matching map...
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::SetTrueGenerationInfo(art::Event const & evt, const int pfpIndex, const bool visibleMode, 
    int &trueGeneration, int &trueParentTrackID, int &trueParentPDG)
{
    trueGeneration = kDefInt;
    trueParentTrackID = kDefInt;
    trueParentPDG = kDefInt;

    // Make sure we have a true neutrino interaction
    if (fNuTrackID == kDefInt)
        return;

    int trueGen = 1; // 1 corresponds to the neutrino
    int currentTrackID = fRecoPFPTrueTrackID[pfpIndex]; // Get the matched track ID
    bool foundParent = false;

    do
    {
        if (fMCParticleMap.find(currentTrackID) == fMCParticleMap.end())
        {
            trueGeneration = kDefInt;
            trueParentTrackID = kDefInt;
            trueParentPDG = kDefInt;
            return;
        }

        const art::Ptr<simb::MCParticle> currentMCParticle = fMCParticleMap.at(currentTrackID);
        currentTrackID = currentMCParticle->Mother();

        // If parent is not reconstructed, move on
        if (visibleMode && (currentTrackID != fNuTrackID))
        {
            if (std::find(fReconstructedMCParticles.begin(), fReconstructedMCParticles.end(), currentTrackID) == 
                fReconstructedMCParticles.end())
            {
                continue;
            }
        }

        if (!foundParent)
        {
            trueParentTrackID = (currentTrackID == fNuTrackID ? fNuTrackID : fMCParticleMap.at(currentTrackID)->TrackId());
            trueParentPDG = (currentTrackID == fNuTrackID ? -1 :fMCParticleMap.at(currentTrackID)->PdgCode());
            foundParent = true;
        }

        ++trueGen;
    }
    while (currentTrackID != fNuTrackID);

    trueGeneration = trueGen;
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoTrackInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpIndex)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fRecoModuleLabel, fTrackModuleLabel))
        return;

    fRecoPFPTrackFitSuccess[pfpIndex] = 1;

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

    if (fRecoTrackRecoEndZ[pfpIndex] > fRecoTrackRecoStartZ[pfpIndex])
    {
        fRecoTrackRecoUpstreamX[pfpIndex] = fRecoTrackRecoStartX[pfpIndex];
        fRecoTrackRecoUpstreamY[pfpIndex] = fRecoTrackRecoStartY[pfpIndex];
        fRecoTrackRecoUpstreamZ[pfpIndex] = fRecoTrackRecoStartZ[pfpIndex];
        fRecoTrackRecoDownstreamX[pfpIndex] = fRecoTrackRecoEndX[pfpIndex];
        fRecoTrackRecoDownstreamY[pfpIndex] = fRecoTrackRecoEndY[pfpIndex];
        fRecoTrackRecoDownstreamZ[pfpIndex] = fRecoTrackRecoEndZ[pfpIndex];
    }
    else
    {
        fRecoTrackRecoDownstreamX[pfpIndex] = fRecoTrackRecoStartX[pfpIndex];
        fRecoTrackRecoDownstreamY[pfpIndex] = fRecoTrackRecoStartY[pfpIndex];
        fRecoTrackRecoDownstreamZ[pfpIndex] = fRecoTrackRecoStartZ[pfpIndex];
        fRecoTrackRecoUpstreamX[pfpIndex] = fRecoTrackRecoEndX[pfpIndex];
        fRecoTrackRecoUpstreamY[pfpIndex] = fRecoTrackRecoEndY[pfpIndex];
        fRecoTrackRecoUpstreamZ[pfpIndex] = fRecoTrackRecoEndZ[pfpIndex];
    }

    fRecoTrackRecoLength[pfpIndex] = track->Length();

    TVector3 upstream_end(fRecoTrackRecoUpstreamX[pfpIndex], fRecoTrackRecoUpstreamY[pfpIndex], fRecoTrackRecoUpstreamZ[pfpIndex]);
    TVector3 downstream_end(fRecoTrackRecoDownstreamX[pfpIndex], fRecoTrackRecoDownstreamY[pfpIndex], fRecoTrackRecoDownstreamZ[pfpIndex]);
    TVector3 vertex_pos(fRecoPFPRecoVertexX[pfpIndex], fRecoPFPRecoVertexY[pfpIndex], fRecoPFPRecoVertexZ[pfpIndex]);

    if ((vertex_pos - upstream_end).Mag() < (vertex_pos - downstream_end).Mag())
    {
        fRecoTrackRecoEndClosestToVertexX[pfpIndex] = fRecoTrackRecoUpstreamX[pfpIndex];
        fRecoTrackRecoEndClosestToVertexY[pfpIndex] = fRecoTrackRecoUpstreamY[pfpIndex];
        fRecoTrackRecoEndClosestToVertexZ[pfpIndex] = fRecoTrackRecoUpstreamZ[pfpIndex];
    }
    else
    {
        fRecoTrackRecoEndClosestToVertexX[pfpIndex] = fRecoTrackRecoDownstreamX[pfpIndex];
        fRecoTrackRecoEndClosestToVertexY[pfpIndex] = fRecoTrackRecoDownstreamY[pfpIndex];
        fRecoTrackRecoEndClosestToVertexZ[pfpIndex] = fRecoTrackRecoDownstreamZ[pfpIndex];
    }

    // Fill momentum variables & energy information
    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(track, evt)));
    fRecoTrackRecoContained[pfpIndex] = energyRecoHandle->longestTrackContained;
    fRecoTrackRecoMomMethod[pfpIndex] = energyRecoHandle->trackMomMethod;

    if (energyRecoHandle->trackMomMethod == 1)
        fRecoTrackRecoMomRange[pfpIndex] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    else if (energyRecoHandle->trackMomMethod == 0)
        fRecoTrackRecoMomMCS[pfpIndex] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    fNumuRecoENu[pfpIndex] = energyRecoHandle->fNuLorentzVector.E();
    fNumuRecoEHad[pfpIndex] = energyRecoHandle->fHadLorentzVector.E();
    fNumuRecoMomLep[pfpIndex] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    // Pandizzle variables
    FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(track, evt));

    /*
    std::cout << "/////////////" << std::endl;
    std::cout << "kMichelNHits: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits) << std::endl; 
    std::cout << "kMichelElectronMVA: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA) << std::endl;
    std::cout << "kMichelRecoEnergyPlane2: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2) << std::endl;
    std::cout << "kTrackDeflecAngleSD: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD) << std::endl;
    std::cout << "kTrackLength: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength) << std::endl;
    std::cout << "kEvalRatio: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio) << std::endl;
    std::cout << "kConcentration: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration) << std::endl;
    std::cout << "kCoreHaloRatio: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio) << std::endl;
    std::cout << "kConicalness: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness) << std::endl;
    std::cout << "kdEdxStart: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart) << std::endl;
    std::cout << "kdEdxEnd: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd) << std::endl;
    std::cout << "kdEdxEndRatio: " << pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio) << std::endl;
    */

    fRecoTrackMichelNHits[pfpIndex] = (double)pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits);
    fRecoTrackMichelElectronMVA[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA);
    fRecoTrackMichelRecoEnergyPlane2[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2);
    fRecoTrackDeflecAngleSD[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
    fRecoTrackLength[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
    fRecoTrackEvalRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
    fRecoTrackConcentration[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
    fRecoTrackCoreHaloRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
    fRecoTrackConicalness[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
    fRecoTrackdEdxStart[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
    fRecoTrackdEdxEnd[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
    fRecoTrackdEdxEndRatio[pfpIndex] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);
    fRecoTrackPandizzleVar[pfpIndex] = pandizzleRecord.GetMVAScore();
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoShowerInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> &pfp,
    const int pfpIndex)
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel))
        return;

    fRecoPFPShowerFitSuccess[pfpIndex] = 1;

    art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, fRecoModuleLabel, fShowerModuleLabel);

    // General shower variables
    fRecoShowerRecoDirX[pfpIndex] = shower->Direction().X();
    fRecoShowerRecoDirY[pfpIndex] = shower->Direction().Y();
    fRecoShowerRecoDirZ[pfpIndex] = shower->Direction().Z();
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
    fRecoShowerRecoMom[pfpIndex] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    if (shower->dEdx().size() > 0)
    {
      for (int i_plane = 0; i_plane < 3; i_plane++)
      {
        fRecoShowerRecodEdx[pfpIndex][i_plane] = shower->dEdx()[i_plane];
        fRecoShowerRecoEnergy[pfpIndex][i_plane] = shower->Energy()[i_plane];
      }
    }

    fNueRecoENu[pfpIndex] = energyRecoHandle->fNuLorentzVector.E();
    fNueRecoEHad[pfpIndex] = energyRecoHandle->fHadLorentzVector.E();
    fNueRecoMomLep[pfpIndex] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    // Pandrizzle
    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower, evt));

    fRecoShowerPandrizzleEvalRatio[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
    fRecoShowerPandrizzleConcentration[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
    fRecoShowerPandrizzleCoreHaloRatio[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
    fRecoShowerPandrizzleConicalness[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
    fRecoShowerPandrizzledEdxBestPlane[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
    fRecoShowerPandrizzleDisplacement[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
    fRecoShowerPandrizzleDCA[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
    fRecoShowerPandrizzleWideness[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
    fRecoShowerPandrizzleEnergyDensity[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);
    fRecoShowerPandrizzlePathwayLengthMin[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
    fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
    fRecoShowerPandrizzleMaxNPostShowerStartHits[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
    fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[pfpIndex] = 
      pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
    fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
    fRecoShowerPandrizzleMaxFoundHitRatio[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
    fRecoShowerPandrizzleMaxInitialGapSize[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
    fRecoShowerPandrizzleMinLargestProjectedGapSize[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
    fRecoShowerPandrizzleNViewsWithAmbiguousHits[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
    fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);
    fRecoShowerPandrizzleModularPathwayLength[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kModularShowerPathwayLengthMin);
    fRecoShowerPandrizzleModularNuVertexChargeWeightedMeanRadialDistance[pfpIndex] =
        pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kModularShowerMaxNuVertexChargeWeightedMeanRadialDistance);
    fRecoShowerPandrizzleModularMaxNShowerHits[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kModularShowerMaxNShowerHits);
    fRecoShowerPandrizzleBDTMethod[pfpIndex] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);

    float pandrizzleScore(pandrizzleRecord.GetMVAScore());
    fRecoShowerBackupPandrizzleScore[pfpIndex] = (std::fabs(fRecoShowerPandrizzleBDTMethod[pfpIndex] - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerEnhancedPandrizzleScore[pfpIndex] = (std::fabs(fRecoShowerPandrizzleBDTMethod[pfpIndex] - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerPandrizzleIsFilled[pfpIndex] = pandrizzleRecord.IsFilled();
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillParentChildLinkInfo(art::Event const & evt)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fRecoModuleLabel))
        return;

    std::vector<art::Ptr<recob::PFParticle>> parentPFPs = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);
    std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fRecoModuleLabel);

    int parentPFPIndex = -1;
    int linkIndex = -1;

    // Loop over all particles as parents
    for (art::Ptr<recob::PFParticle> parentPFP : parentPFPs)
    {
        if ((std::fabs(parentPFP->PdgCode()) == 12) || (std::fabs(parentPFP->PdgCode()) == 14) || (std::fabs(parentPFP->PdgCode()) == 16))
            continue;

        parentPFPIndex++;

        if (parentPFPIndex >= kMaxPFParticles)
            break;

        // Loop over all particles as children
        int childPFPIndex = -1;

        for (art::Ptr<recob::PFParticle> childPFP : childPFPs)
        {
            if ((std::fabs(childPFP->PdgCode()) == 12) || (std::fabs(childPFP->PdgCode()) == 14) || (std::fabs(childPFP->PdgCode()) == 16))
                continue;

            childPFPIndex++;

            if (childPFPIndex >= kMaxPFParticles)
                break;

            // Increase number of links
            linkIndex++;
            fNParentChildLinks++;

            // Index information   
            fParentPFPIndex[linkIndex] = parentPFPIndex;
            fChildPFPIndex[linkIndex] = childPFPIndex;

            // Fill true parent-child link info
            FillTrueParentChildLinkInfo(linkIndex, parentPFPIndex, childPFPIndex);

            // Fill reco parent-child link info
            FillRecoParentChildLinkInfo(evt, childPFP, parentPFP, linkIndex);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillTrueParentChildLinkInfo(const int linkIndex, const int parentPFPIndex, 
    const int childPFPIndex)
{
    if (fRecoPFPTrueTrackID[parentPFPIndex] < 0) // no reco match for parent
    {
        fTrueParentChildLink[linkIndex] = false;
    }
    else if (fRecoPFPTrueVisibleParentTrackID[childPFPIndex] < 0) // no true parent identified 
    {
        fTrueParentChildLink[linkIndex] = false;
    }
    else if (fRecoPFPTrueVisibleParentTrackID[childPFPIndex] == fRecoPFPTrueTrackID[parentPFPIndex]) // correct child link
    {
        fTrueParentChildLink[linkIndex] = true;
    }
    else
    {
        fTrueParentChildLink[linkIndex] = false;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillRecoParentChildLinkInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> childPFP, 
    art::Ptr<recob::PFParticle> parentPFP, const int linkIndex)
{ 
    // Parent information
    fParentTrackScore[linkIndex] = HierarchyUtils::GetTrackScore(evt, parentPFP, fRecoModuleLabel);
    fParentBraggVariable[linkIndex] = HierarchyUtils::GetBraggVariable();
    fParentEndRegionNHits[linkIndex] = HierarchyUtils::GetEndRegionNHits();
    fParentEndRegionNParticles[linkIndex] = HierarchyUtils::GetEndRegionNParticles();
    fParentEndRegionRToWall[linkIndex] = HierarchyUtils::GetEndRegionRToWall();

    // Edge information
    fVertexSeparation[linkIndex] = HierarchyUtils::GetVertexSeparation();
    fSeparation3D[linkIndex] = HierarchyUtils::GetSeparation3D();
    fEnergyRatio[linkIndex] = HierarchyUtils::GetEnergyRatio();
    fPIDLinkType[linkIndex] = HierarchyUtils::GetPIDLinkType();
    fOpeningAngle[linkIndex] = HierarchyUtils::GetOpeningAngle();
    fTrackShowerLinkType[linkIndex] = HierarchyUtils::GetTrackShowerLinkType();
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunTrackSelection(art::Event const & evt)
{
    // Pandizzle
    RunPandizzleTrackSelection();

    // DeepPan
    RunDeepPanTrackSelection();

    // Ivysaurus
    RunIvysaurusTrackSelection();

    // Longest Length
    RunLongestLengthTrackSelection();
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunPandizzleTrackSelection()
{
    float highestPandizzleScore = -1.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoTrackPandizzleVar[i] < highestPandizzleScore)
            continue;

        highestPandizzleScore = fRecoTrackPandizzleVar[i];
        fSelTrackPandizzleSelf = fRecoPFPSelf[i];
        fSelTrackPandizzleIndex = i;
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunDeepPanTrackSelection()
{
    float highestDeepPanScore = -0.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoPFPDeepPanMuVar[i] < highestDeepPanScore)
            continue;

        highestDeepPanScore = fRecoPFPDeepPanMuVar[i];
        fSelTrackDeepPanSelf = fRecoPFPSelf[i];
        fSelTrackDeepPanIndex = i;
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunIvysaurusTrackSelection()
{
    float highestIvysaurusScore = -0.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoPFPIvysaurusParticleType[i] != 0)
            continue;

        if (fRecoPFPIvysaurusMuon[i] < highestIvysaurusScore)
            continue;

        highestIvysaurusScore = fRecoPFPIvysaurusMuon[i];
        fSelTrackIvysaurusSelf = fRecoPFPSelf[i];
        fSelTrackIvysaurusIndex = i;
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunLongestLengthTrackSelection()
{
    float longestLength = -0.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoTrackRecoLength[i] < longestLength)
            continue;

        longestLength = fRecoTrackRecoLength[i];
        fSelTrackLongestLengthSelf = fRecoPFPSelf[i];
        fSelTrackLongestLengthIndex = i;
    }
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunShowerSelection(art::Event const & evt)
{
    // Pandrizzle
    RunPandrizzleShowerSelection();

    // Ivysaurus
    RunIvysaurusShowerSelection();

    // HighestEnergy
    RunHighestEnergyShowerSelection();
}

//////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunPandrizzleShowerSelection()
{
    bool foundEnhanced = false;
    float highestEnhancedPandrizzleScore = -1.1;
    float highestBackupPandrizzleScore = -1.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if ((std::fabs(fRecoShowerPandrizzleBDTMethod[i] - 2.0) < std::numeric_limits<float>::epsilon()) && 
            (fRecoShowerEnhancedPandrizzleScore[i] > highestEnhancedPandrizzleScore))
        {
            foundEnhanced = true;

            highestEnhancedPandrizzleScore = fRecoShowerEnhancedPandrizzleScore[i];
            fSelShowerPandrizzleSelf = fRecoPFPSelf[i];
            fSelShowerPandrizzleIndex = i;
        }

        if (foundEnhanced)
            continue;

        if ((std::fabs(fRecoShowerPandrizzleBDTMethod[i] - 1.0) < std::numeric_limits<float>::epsilon()) && 
            (fRecoShowerBackupPandrizzleScore[i] > highestBackupPandrizzleScore))
        {
            highestBackupPandrizzleScore = fRecoShowerBackupPandrizzleScore[i];
            fSelShowerPandrizzleSelf = fRecoPFPSelf[i];
            fSelShowerPandrizzleIndex = i;
        }
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunIvysaurusShowerSelection()
{
    float highestIvysaurusScore = -0.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoPFPIvysaurusParticleType[i] != 3)
            continue;

        if (fRecoPFPIvysaurusElectron[i] < highestIvysaurusScore)
            continue;

        highestIvysaurusScore = fRecoPFPIvysaurusElectron[i];
        fSelShowerIvysaurusSelf = fRecoPFPSelf[i];
        fSelShowerIvysaurusIndex = i;
    }
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunHighestEnergyShowerSelection()
{
    float highestEnergy = -0.1;

    for (int i = 0; i < fNRecoPFPs; ++i)
    {
        if (!fRecoPFPIsPrimary[i])
            continue;

        if (fRecoPFPSelf[i] < 0)
            continue;

        if (fRecoShowerRecoEnergy[i][2] < highestEnergy)
            continue;

        highestEnergy = fRecoShowerRecoEnergy[i][2];
        fSelShowerHighestEnergySelf = fRecoPFPSelf[i];
        fSelShowerHighestEnergyIndex = i;
    }
}

//////////////////////////////////////////////////////////////////

TVector3 FDSelection::CCNuSelection::ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector){
  TVector3 projected_vector = vector_to_project - (vector_to_project.Dot(plane_norm_vector) / (plane_norm_vector.Mag() * plane_norm_vector.Mag()))*plane_norm_vector;
  return projected_vector;
}

//////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(FDSelection::CCNuSelection)
