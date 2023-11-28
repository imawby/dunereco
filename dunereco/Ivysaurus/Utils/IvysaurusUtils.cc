#include "IvysaurusUtils.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

namespace IvysaurusUtils
{

/////////////////////////////////////////////////////////////

const TVector3 ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
    const PandoraView hitType)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const geo::WireID hitWireID = hit->WireID();
    const geo::CryostatID cryostatID(hitWireID.Cryostat);
    const double hitTime = hit->PeakTime();
    const double xCoord = detProp.ConvertTicksToX(hitTime, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);

    // Get hit Y and Z coordinates, based on central position of wire
    geo::Point_t xyz = theGeometry->Cryostat(cryostatID).TPC(hitWireID.TPC).Plane(hitWireID.Plane).Wire(hitWireID.Wire).GetCenter();

    return TVector3(xCoord, 0.f, hitType == TPC_VIEW_U ? YZToU(xyz.Y(), xyz.Z()) : hitType == TPC_VIEW_V ? YZToV(xyz.Y(), xyz.Z()) : YZToW(xyz.Y(), xyz.Z()));
}

/////////////////////////////////////////////////////////////

const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const PandoraView pandoraView)
{
    const float xCoord = inputPosition3D.X();
    const float yCoord = inputPosition3D.Y();
    const float zCoord = inputPosition3D.Z();

    return TVector3(xCoord, 0.f, pandoraView == TPC_VIEW_U ? YZToU(yCoord, zCoord) : pandoraView == TPC_VIEW_V ? YZToV(yCoord, zCoord) : YZToW(yCoord, zCoord));
}


/////////////////////////////////////////////////////////////

const PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit)
{
    const geo::WireID hitWireID(hit->WireID());
    const geo::View_t hitView(hit->View());
    const geo::View_t thisPandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hitWireID.Cryostat, hitWireID.TPC, hitView));

    if (thisPandoraView == geo::kW || thisPandoraView == geo::kY)
        return TPC_VIEW_W;
    else if (thisPandoraView == geo::kU)
        return TPC_VIEW_U;
    else if (thisPandoraView == geo::kV)
        return TPC_VIEW_V;
    else
        throw cet::exception("ivysaur::GridManager") << "wire view not recognised";
}

/////////////////////////////////////////////////////////////

float YZToU(const float yCoord, const float zCoord)
{
    const float m_uWireAngle = 0.623257100582;

    return (zCoord * std::cos(m_uWireAngle)) - (yCoord * std::sin(m_uWireAngle));
}

/////////////////////////////////////////////////////////////

float YZToV(const float yCoord, const float zCoord)
{
    const float m_vWireAngle = -0.623257100582;

    return (zCoord * std::cos(m_vWireAngle)) - (yCoord * std::sin(m_vWireAngle));
}

/////////////////////////////////////////////////////////////

float YZToW(const float yCoord, const float zCoord)
{
    const float m_wWireAngle = 0.0;

    return (zCoord * std::cos(m_wWireAngle)) - (yCoord * std::sin(m_wWireAngle));
}

/////////////////////////////////////////////////////////////

double CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits, 
  const std::vector<art::Ptr<recob::Hit>> &eventHits, const int trackID)
{
    int nMatchesInSelHits = 0;
    int nMatchesInAllHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1); 
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);;
        if (matchedID == trackID) nMatchesInAllHits++;
    }

    const double completeness = (nMatchesInAllHits > 0) ? static_cast<double>(nMatchesInSelHits) / static_cast<double>(nMatchesInAllHits) : 0.0;

    return completeness;
}

/////////////////////////////////////////////////////////////////////////

double HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits,
    const int trackID)
{
    int nMatchesInSelHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    const double purity = (selectedHits.size() > 0) ? static_cast<double>(nMatchesInSelHits) / static_cast<double>(selectedHits.size()) : 0.0;

    std::cout << "purity: " << purity << std::endl;

    return purity;
}
}




