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

#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

class GridManager
{
  public:
    enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

  class Grid
  {
    public:
      Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, const PandoraView pandoraView, const bool isInitialised);

      unsigned int GetAxisDimensions() const;
      std::vector<float> GetDriftBoundaries() const;
      std::vector<float> GetWireBoundaries() const;
      std::vector<std::vector<float>> GetGridValues() const;
      std::vector<std::vector<float>> GetCountValues() const;
      PandoraView GetPandoraView() const;
      bool IsInitialised() const;
      bool IsNormalised() const;
      bool IsInsideGrid(const TVector3 &position) const;
      void AddToGrid(const TVector3 &position, const float energy, const float weight);
      void NormaliseGrid();

    private:
      unsigned int m_axisDimensions; // number of bins on each axis
      std::vector<float> m_driftBoundaries;
      std::vector<float> m_wireBoundaries;
      std::vector<std::vector<float>> m_gridValues; // driftBin : [wireBins]
      std::vector<std::vector<float>> m_countValues; // driftBin : [wireBins]
      GridManager::PandoraView m_pandoraView;
      bool m_isInitialised;
      bool m_isNormalised;
  };

  public:
    GridManager(const fhicl::ParameterSet& pset);
    ~GridManager();

    // Function to place the grid in space
    Grid ObtainViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const PandoraView tpcView) const;

    // Function to fill the grid
    void FillViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        GridManager::Grid &grid) const;

    // Function to obtain the Pandora view of a LArSoft hit
    const GridManager::PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit) const;

    // Function to find the 2D Pandora coordinate of a LArSoft hit
    const TVector3 ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
        const GridManager::PandoraView hitType) const;

  private:
    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    // Function to project a 3D coordinate into a specified Pandora 2D view
    const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const GridManager::PandoraView pandoraView) const;

    // Function to obtain the Pandora U coordinate from LArSoft Y/Z coordinates
    float YZToU(const float yCoord, const float zCoord) const;

    // Function to obtain the Pandora V coordinate from LArSoft Y/Z coordinates
    float YZToV(const float yCoord, const float zCoord) const;

    // Function to obtain the Pandora W coordinate from LArSoft Y/Z coordinates
    float YZToW(const float yCoord, const float zCoord) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;

    float m_gridSize3D;
    unsigned int m_dimensions;
    float m_uWireAngle;
    float m_vWireAngle;
    float m_wWireAngle;
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

inline GridManager::PandoraView GridManager::Grid::GetPandoraView() const 
{ 
    return m_pandoraView; 
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsInitialised() const 
{ 
    return m_isInitialised;
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsNormalised() const 
{ 
    return m_isNormalised;
}

}

#endif  // GRIDMANAGER_H
