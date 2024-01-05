////////////////////////////////////////////////////////////////////////
/// \file    GridManager.cxx
/// \brief   A class to manage the Ivysaurus 2D grid input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef GRIDMANAGER_H
#define GRIDMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

class GridManager
{
  public:
  class Grid
  {
    public:
      Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, const float maxGridEntry, 
          const IvysaurusUtils::PandoraView pandoraView, const bool isInitialised);

      unsigned int GetAxisDimensions() const;
      std::vector<float> GetDriftBoundaries() const;
      std::vector<float> GetWireBoundaries() const;
      std::vector<std::vector<float>> GetGridValues() const;
      std::vector<std::vector<float>> GetCountValues() const;
      IvysaurusUtils::PandoraView GetPandoraView() const;
      bool IsInitialised() const;
      bool IsAveraged() const;
      bool IsNormalised() const;
      bool IsInsideGrid(const TVector3 &position) const;
      void AddToGrid(const TVector3 &position, const float energy, const float weight);
      void AverageGrid();
      void NormaliseGrid();

    private:
      unsigned int m_axisDimensions; // number of bins on each axis
      float m_maxGridEntry;
      std::vector<float> m_driftBoundaries;
      std::vector<float> m_wireBoundaries;
      std::vector<std::vector<float>> m_gridValues; // driftBin : [wireBins]
      std::vector<std::vector<float>> m_countValues; // driftBin : [wireBins]
      IvysaurusUtils::PandoraView m_pandoraView;
      bool m_isInitialised;
      bool m_isAveraged;
      bool m_isNormalised;
  };

    GridManager(const fhicl::ParameterSet& pset);
    ~GridManager();

    // Function to place the grid in space
    Grid ObtainViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const IvysaurusUtils::PandoraView pandoraView, const bool isStart) const;

    // Function to fill the grid
    void FillViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        GridManager::Grid &grid) const;

  private:
    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    bool GetStartExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    bool GetStartExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    std::string m_showerModuleLabel;

    float m_gridSize3D;
    unsigned int m_dimensions;
    float m_maxGridEntry;
    float m_recombFactor;
    calo::CalorimetryAlg m_calorimetryAlg;
};

/////////////////////////////////////////////////////////////

inline unsigned int GridManager::Grid::GetAxisDimensions() const 
{ 
    return m_axisDimensions; 
}

/////////////////////////////////////////////////////////////

inline std::vector<float> GridManager::Grid::GetDriftBoundaries() const 
{ 
    return m_driftBoundaries; 
}

/////////////////////////////////////////////////////////////

inline std::vector<float> GridManager::Grid::GetWireBoundaries() const 
{ 
    return m_wireBoundaries; 
}

/////////////////////////////////////////////////////////////

inline std::vector<std::vector<float>> GridManager::Grid::GetGridValues() const 
{ 
    return m_gridValues; 
}

/////////////////////////////////////////////////////////////

inline std::vector<std::vector<float>> GridManager::Grid::GetCountValues() const 
{ 
    return m_countValues; 
}

/////////////////////////////////////////////////////////////

inline IvysaurusUtils::PandoraView GridManager::Grid::GetPandoraView() const 
{ 
    return m_pandoraView; 
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsInitialised() const 
{ 
    return m_isInitialised;
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsAveraged() const 
{ 
    return m_isAveraged;
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsNormalised() const 
{ 
    return m_isNormalised;
}

}

#endif  // GRIDMANAGER_H
