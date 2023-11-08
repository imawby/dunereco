////////////////////////////////////////////////////////////////////////
/// \file    TrackVarManager.cxx
/// \brief   A class to manage the Ivysaurus 2D track variable input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/Ivysaurus/Managers/TrackVarManager.h"


namespace ivysaurus
{

TrackVarManager::TrackVars::TrackVars() : 
        m_nTrackChildren(-1),
        m_nShowerChildren(-1),
        m_nGrandChildren(-1),
        m_trackLength(-1.f),
        m_wobble(-1.f)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TrackVarManager::TrackVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel"))
{
}

/////////////////////////////////////////////////////////////

TrackVarManager::~TrackVarManager()
{
}

/////////////////////////////////////////////////////////////

bool TrackVarManager::EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TrackVarManager::TrackVars &trackVars) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    FillHierarchyInfo(evt, pfparticle, trackVars);
    FillTrackLength(track, trackVars);
    FillWobble(track, trackVars);

    return true;
}

/////////////////////////////////////////////////////////////

void TrackVarManager::FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TrackVarManager::TrackVars &trackVars) const
{
    const std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);

    int nTracks = 0;
    int nShowers = 0;
    int nGrandChildren = 0;

    for (const art::Ptr<recob::PFParticle> &childPFP : childPFPs)
    {
        if (childPFP->PdgCode() == 13)
            ++nTracks;
        else if (childPFP->PdgCode() == 11)
            ++nShowers;

        const std::vector<art::Ptr<recob::PFParticle>> grandChildPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(childPFP, evt, m_recoModuleLabel);

        nGrandChildren += grandChildPFPs.size();
    }

    trackVars.SetNTrackChildren(nTracks);
    trackVars.SetNShowerChildren(nShowers);
    trackVars.SetNGrandChildren(nGrandChildren);
}

/////////////////////////////////////////////////////////////

void TrackVarManager::FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const
{
    const float trackLength = track->Length();

    trackVars.SetTrackLength(trackLength);
}

/////////////////////////////////////////////////////////////
// This is copied from Dom's Pandizzle module
void TrackVarManager::FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const
{
  // This follows the MicroBooNE method for MCS

  //Get the number of points
  size_t NPoints = track->NumberTrajectoryPoints();
  //Store the directions between adjacent points on a vector
  std::vector<TVector3> directions;
  for (size_t i_point = 0; i_point < NPoints-1; i_point++){
    TVector3 position_i(track->TrajectoryPoint(i_point).position.X(), track->TrajectoryPoint(i_point).position.Y(), track->TrajectoryPoint(i_point).position.Z());
    TVector3 position_iplus1(track->TrajectoryPoint(i_point+1).position.X(), track->TrajectoryPoint(i_point+1).position.Y(), track->TrajectoryPoint(i_point+1).position.Z());
    TVector3 direction = (position_iplus1-position_i).Unit();
    directions.push_back(direction);
  }

  //Loop through the direction and compare adjacent elements
  std::vector<double> deflection_angles;
  for (size_t i_dir = 0; i_dir < directions.size()-1; i_dir++){
    //Aim: rotate both direction so that the first direction is parallel to the z-axis.  
    //Then take the x-projection of scattered track and calculate the angle between that and the first direction
    TVector3 z_axis(0,0,1);
    TVector3 direction_first = directions[i_dir];
    TVector3 direction_second = directions[i_dir+1];

    //Ignore if either direction is 0 (i.e. not 1)
    if (direction_first.Mag() < 0.999 || direction_second.Mag() < 0.999){
      continue;
    }
    double angle_dir_first_z_axis = direction_first.Angle(z_axis);
    TVector3 orthogonal_vector = direction_first.Cross(z_axis);
    if (orthogonal_vector.Unit().Mag() < 0.999){
        continue;
    }
    direction_first.Rotate(angle_dir_first_z_axis, orthogonal_vector);
    direction_second.Rotate(angle_dir_first_z_axis, orthogonal_vector);

    //Now work out the angle between the vectors in the x-z plane
    direction_first.SetY(0);
    direction_second.SetY(0);
    double dot_product = direction_first.Dot(direction_second);
    dot_product = std::min(std::max(dot_product,-1.),1.);
    double angle = acos(dot_product) * 180/3.142;

    //define +x as a +angle
    if (direction_second.X() < 0) angle*=-1;
    deflection_angles.push_back(angle);
  }

  double angle_mean = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_mean += deflection_angles[i_angle];
  }

  if (deflection_angles.size()>0) angle_mean/=deflection_angles.size();
  else angle_mean=-100;

  double angle_var = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_var = (deflection_angles[i_angle] - angle_mean)*(deflection_angles[i_angle] - angle_mean);
  }

  if (deflection_angles.size() > 1) angle_var /= (deflection_angles.size()-1);
  else angle_var = -2.;

  if (angle_var > 0.0)
      trackVars.SetWobble(sqrt(angle_var));
}

/////////////////////////////////////////////////////////////

void TrackVarManager::Reset(TrackVarManager::TrackVars &trackVars) const
{
    trackVars.SetNTrackChildren(-1);
    trackVars.SetNShowerChildren(-1);
    trackVars.SetNGrandChildren(-1);
    trackVars.SetTrackLength(-1.f);
    trackVars.SetWobble(-1.f);
}

/////////////////////////////////////////////////////////////

}
