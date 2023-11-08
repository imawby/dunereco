////////////////////////////////////////////////////////////////////////
/// \file    ShowerVarManager.h
/// \brief   A class to manage the Ivysaurus shower variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef SHOWERMANAGER_H
#define SHOWERMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

namespace ivysaurus
{

class ShowerVarManager
{
  public:
  class ShowerVars
  {
    public:
      ShowerVars();

      float GetDisplacement() const;
      float GetInitialGapSize() const;
      float GetLargestGapSize() const;
      float GetPathwayLength() const;
      float GetPathwayScatteringAngle2D() const;
      float GetNShowerHits() const;
      float GetFoundHitRatio() const;
      float GetScatterAngle() const;
      float GetOpeningAngle() const;
      float GetNuVertexEnergyAsymmetry() const;
      float GetNuVertexEnergyWeightedMeanRadialDistance() const;
      float GetShowerStartEnergyAsymmetry() const;
      float GetShowerStartMoliereRadius() const;
      float GetNAmbiguousViews() const;
      float GetUnaccountedEnergy() const;
      ////
      void SetDisplacement(const float displacement);
      void SetInitialGapSize(const float initialGapSize);
      void SetLargestGapSize(const float largestGapSize);
      void SetPathwayLength(const float pathwayLength);
      void SetPathwayScatteringAngle2D(const float pathwayScatteringAngle2D);
      void SetNShowerHits(const float nShowerHits);
      void SetFoundHitRatio(const float foundHitRatio);
      void SetScatterAngle(const float scatterAngle);
      void SetOpeningAngle(const float openingAngle);
      void SetNuVertexEnergyAsymmetry(const float nuVertexEnergyAsymmetry);
      void SetNuVertexEnergyWeightedMeanRadialDistance(const float nuVertexEnergyWeightedMeanRadialDistance);
      void SetShowerStartEnergyAsymmetry(const float showerStartEnergyAsymmetry);
      void SetShowerStartMoliereRadius(const float showerStartMoliereRadius);
      void SetNAmbiguousViews(const float nAmbiguousViews);
      void SetUnaccountedEnergy(const float unaccountedEnergy);

    private:
      float m_displacement;
      float m_initialGapSize;
      float m_largestGapSize;
      float m_pathwayLength;
      float m_pathwayScatteringAngle2D;
      float m_nShowerHits;
      float m_foundHitRatio;
      float m_scatterAngle;
      float m_openingAngle;
      float m_nuVertexEnergyAsymmetry;
      float m_nuVertexEnergyWeightedMeanRadialDistance;
      float m_showerStartEnergyAsymmetry;
      float m_showerStartMoliereRadius;
      float m_nAmbiguousViews;
      float m_unaccountedEnergy;
  };

    ShowerVarManager(const fhicl::ParameterSet& pset);
    ~ShowerVarManager();

    bool EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, ShowerVarManager::ShowerVars &showerVars) const;
    void Reset(ShowerVarManager::ShowerVars &showerVars) const;

  private:
    void FillDisplacement(const art::Event &evt, const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;

    void FillConnectionPathwayVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        ShowerVarManager::ShowerVars &showerVars) const;

    std::string m_recoModuleLabel;
    std::string m_showerModuleLabel;
};

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetDisplacement() const
{
    return m_displacement;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetInitialGapSize() const
{
    return m_initialGapSize;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetLargestGapSize() const
{
    return m_largestGapSize;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetPathwayLength() const
{
    return m_pathwayLength;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetPathwayScatteringAngle2D() const
{
    return m_pathwayScatteringAngle2D;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNShowerHits() const
{
    return m_nShowerHits;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetFoundHitRatio() const
{
    return m_foundHitRatio;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetScatterAngle() const
{
    return m_scatterAngle;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetOpeningAngle() const
{
    return m_openingAngle;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexEnergyAsymmetry() const
{
    return m_nuVertexEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexEnergyWeightedMeanRadialDistance() const
{
    return m_nuVertexEnergyWeightedMeanRadialDistance;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetShowerStartEnergyAsymmetry() const
{
    return m_showerStartEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetShowerStartMoliereRadius() const
{
    return m_showerStartMoliereRadius;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNAmbiguousViews() const
{
    return m_nAmbiguousViews;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetUnaccountedEnergy() const
{
    return m_unaccountedEnergy;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetDisplacement(const float displacement) 
{
    m_displacement = displacement;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetInitialGapSize(const float initialGapSize) 
{
    m_initialGapSize = initialGapSize;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetLargestGapSize(const float largestGapSize) 
{
    m_largestGapSize = largestGapSize;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetPathwayLength(const float pathwayLength) 
{
    m_pathwayLength = pathwayLength;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetPathwayScatteringAngle2D(const float pathwayScatteringAngle2D) 
{
    m_pathwayScatteringAngle2D = pathwayScatteringAngle2D;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNShowerHits(const float nShowerHits) 
{
    m_nShowerHits = nShowerHits;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetFoundHitRatio(const float foundHitRatio) 
{
    m_foundHitRatio =  foundHitRatio;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetScatterAngle(const float scatterAngle) 
{
    m_scatterAngle = scatterAngle;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetOpeningAngle(const float openingAngle) 
{
    m_openingAngle = openingAngle;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexEnergyAsymmetry(const float nuVertexEnergyAsymmetry) 
{
    m_nuVertexEnergyAsymmetry = nuVertexEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexEnergyWeightedMeanRadialDistance(const float nuVertexEnergyWeightedMeanRadialDistance) 
{
    m_nuVertexEnergyWeightedMeanRadialDistance = nuVertexEnergyWeightedMeanRadialDistance;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetShowerStartEnergyAsymmetry(const float showerStartEnergyAsymmetry) 
{
    m_showerStartEnergyAsymmetry = showerStartEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetShowerStartMoliereRadius(const float showerStartMoliereRadius) 
{
    m_showerStartMoliereRadius = showerStartMoliereRadius;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNAmbiguousViews(const float nAmbiguousViews) 
{
    m_nAmbiguousViews = nAmbiguousViews;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetUnaccountedEnergy(const float unaccountedEnergy) 
{
    m_unaccountedEnergy = unaccountedEnergy;
}

/////////////////////////////////////////////////////////////



}

#endif  // SHOWERVARMANAGER_H
