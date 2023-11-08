////////////////////////////////////////////////////////////////////////
/// \file    TrackVarManager.h
/// \brief   A class to manage the Ivysaurus track variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef TRACKMANAGER_H
#define TRACKMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

namespace ivysaurus
{

class TrackVarManager
{
  public:
  class TrackVars
  {
    public:
      TrackVars();

      int GetNTrackChildren() const;
      int GetNShowerChildren() const;
      int GetNGrandChildren() const;
      float GetTrackLength() const;
      float GetWobble() const;

      void SetNTrackChildren(const int nTrackChildren);
      void SetNShowerChildren(const int nShowerChildren);
      void SetNGrandChildren(const int nGrandChildren);
      void SetTrackLength(const float trackLength);
      void SetWobble(const float wobble);

    private:
      int m_nTrackChildren;
      int m_nShowerChildren;
      int m_nGrandChildren;
      float m_trackLength;
      float m_wobble;
  };

    TrackVarManager(const fhicl::ParameterSet& pset);
    ~TrackVarManager();

    bool EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    void Reset(TrackVarManager::TrackVars &trackVars) const;

  private:
    void FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    void FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    void FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;

    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
};

/////////////////////////////////////////////////////////////

inline int TrackVarManager::TrackVars::GetNTrackChildren() const
{
    return m_nTrackChildren;
}

/////////////////////////////////////////////////////////////

inline int TrackVarManager::TrackVars::GetNShowerChildren() const
{
    return m_nShowerChildren;
}

/////////////////////////////////////////////////////////////

inline int TrackVarManager::TrackVars::GetNGrandChildren() const
{
    return m_nGrandChildren;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetTrackLength() const
{
    return m_trackLength;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetWobble() const
{
    return m_wobble;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNTrackChildren(const int nTrackChildren)
{
    m_nTrackChildren = nTrackChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNShowerChildren(const int nShowerChildren)
{
    m_nShowerChildren = nShowerChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNGrandChildren(const int nGrandChildren)
{
    m_nGrandChildren = nGrandChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetTrackLength(const float trackLength)
{
    m_trackLength = trackLength;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetWobble(const float wobble)
{
    m_wobble = wobble;
}

/////////////////////////////////////////////////////////////




}

#endif  // TRACKVARMANAGER_H
