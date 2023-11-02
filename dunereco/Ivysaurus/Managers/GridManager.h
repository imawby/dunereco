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

namespace ivysaurus
{

class GridManager
{
  enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

  class Grid
  {
    public:
      Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, const PandoraView pandoraView, const bool isInitialised);

      const PandoraView GetPandoraView() const;
      const bool IsInitialised() const;
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
    const Grid ObtainViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const PandoraView tpcView) const;

    // Function to fill the grid
    void FillViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        GridManager::Grid &grid) const;

  private:
    // Function to find the 2D Pandora coordinate of a LArSoft hit
    const TVector3 ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
        const GridManager::PandoraView hitType) const;

    // Function to project a 3D coordinate into a specified Pandora 2D view
    const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const GridManager::PandoraView pandoraView) const;

    // Function to obtain the Pandora U coordinate from LArSoft Y/Z coordinates
    float YZToU(const float yCoord, const float zCoord) const;

    // Function to obtain the Pandora V coordinate from LArSoft Y/Z coordinates
    float YZToV(const float yCoord, const float zCoord) const;

    // Function to obtain the Pandora W coordinate from LArSoft Y/Z coordinates
    float YZToW(const float yCoord, const float zCoord) const;

    // Function to obtain the Pandora view of a LArSoft hit
    const GridManager::PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;

    float m_gridSize3D;
    unsigned int m_dimensions;
    float m_uWireAngle;
    float m_vWireAngle;
    float m_wWireAngle;
};

/////////////////////////////////////////////////////////////

inline const GridManager::PandoraView GridManager::Grid::GetPandoraView() const 
{ 
    return m_pandoraView; 
}

/////////////////////////////////////////////////////////////

inline const bool GridManager::Grid::IsInitialised() const 
{ 
    return m_isInitialised;
}






}

#endif  // GRIDMANAGER_H
