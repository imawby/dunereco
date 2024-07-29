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

namespace HierarchyUtils
{

 double DEFAULT_DOUBLE = -9999.0;
 int DEFAULT_INT = -999;

 // PFP-level
 double GetTrackScore(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel);
 double GetBraggVariable();
 int GetEndRegionNHits();
 int GetEndRegionNParticles();
 double GetEndRegionRToWall();
 // Edge-level
 double GetVertexSeparation();
 double GetSeparation3D();
 double GetEnergyRatio();
 int GetPIDLinkType();
 double GetOpeningAngle();
 int GetTrackShowerLinkType();
}

#endif
