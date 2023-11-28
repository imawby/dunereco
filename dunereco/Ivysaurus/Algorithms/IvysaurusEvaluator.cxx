////////////////////////////////////////////////////////////////////////
/// \file    IvysaurusEvaluator.cxx
/// \brief   Functions to use the Ivysaurus network
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dunereco/Ivysaurus/Algorithms/IvysaurusEvaluator.h"

#include "cetlib/getenv.h"

namespace ivysaurus
{

/////////////////////////////////////////////////////////////

IvysaurusEvaluator::IvysaurusEvaluator(const fhicl::ParameterSet& pset) :
    m_networkDir(pset.get<std::string>("NetworkPath")),
    m_networkName(pset.get<std::string>("NetworkName"))
{
}

/////////////////////////////////////////////////////////////

void IvysaurusEvaluator::IvysaurusUseEvaluate(const art::Ptr<recob::PFParticle> pfparticle, const art::Event &evt) const
{

}

}
