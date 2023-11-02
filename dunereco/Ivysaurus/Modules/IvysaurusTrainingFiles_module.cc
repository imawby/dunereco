/**
 *  @file   dunereco/Ivysaurus/Modules/IvysaurusTrainingFiles_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "dunereco/Ivysaurus/Managers/GridManager.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ivysaurus
{

/**
 *  @brief  IvysaurusTrainingFiles class
 */
class IvysaurusTrainingFiles : public art::EDAnalyzer
{
public:
   IvysaurusTrainingFiles(fhicl::ParameterSet const &pset);
   virtual ~IvysaurusTrainingFiles();

   void beginJob();
   void endJob();
   void analyze(const art::Event &evt);

private:

  GridManager m_gridManager;
  std::string m_recoModuleLabel;
};

DEFINE_ART_MODULE(IvysaurusTrainingFiles)

} // namespace ivysaurus

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <random>

namespace ivysaurus
{

IvysaurusTrainingFiles::IvysaurusTrainingFiles(fhicl::ParameterSet const &pset) : 
  art::EDAnalyzer(pset),
  m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
  m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

IvysaurusTrainingFiles::~IvysaurusTrainingFiles()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::analyze(const art::Event &evt)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

} //namespace ivysaurus

