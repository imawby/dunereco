#ifndef NEWHIERARCHYUTILS_H_SEEN
#define NEWHIERARCHYUTILS_H_SEEN

///////////////////////////////////////////////
// HierarchyUtils.h
//
// Isobel Mawby i.mawby1@lancaster.ac.uk
///////////////////////////////////////////////

// Framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataobj/RecoBase/PFParticle.h"

// c++
#include <vector>
#include <map>

// add private??

namespace HierarchyUtils
{

 double DEFAULT_DOUBLE = -9999.0;
 int DEFAULT_INT = -999;

 // PFP-level
 double GetTrackScore(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel);

 double GetNuVertexSeparation(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, 
     const std::string recoModuleLabel);

 double GetBraggVariable();

 int GetEndRegionNHits(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string trackModuleLabel,  
     const std::string recoModuleLabel, const double separationThreshold);

 int GetEndRegionNParticles(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string trackModuleLabel,  
     const std::string recoModuleLabel, const double separationThreshold);

 double GetEndRegionRToWall(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel, 
     const std::string trackModuleLabel);

 // Edge-level
 double GetVertexSeparation(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
     const std::string recoModuleLabel);

 double GetSeparation3D(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
     const std::string recoModuleLabel);

 void GetLinkConnectionInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
     const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &connectionVars);

 bool IsPandoraApprovedTrack(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel, 
     const std::string trackModuleLabel);

 bool IsInBoundingBox(const TVector3 &boundary1, const TVector3 &boundary2, const TVector3 &testPoint, const float buffer);

 void GetParentConnectionPointVars(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP,
     const std::string recoModuleLabel, const double searchRegion, std::map<std::string, double> &connectionVars);

 void RunPCA(const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints, std::vector<double> &eVals, std::vector<TVector3> &eVecs);

 double GetChargeRatio(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
     const std::string recoModuleLabel);

 double GetPFPCharge(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel);

 int GetPIDLinkTypeWithIvysaurus(const int parentParticleType, const int childParticleType);

 int GetPIDLinkTypeWithPDG(const int parentPDG, const int childPDG);

 double GetOpeningAngle(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel);

 bool GetParticleDirection(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfp, const std::string recoModuleLabel, 
    const float searchRegion, TVector3 &pfpDirection);

 int GetTrackShowerLinkType(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel);

}

#endif
