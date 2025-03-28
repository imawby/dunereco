include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(
    LIB_LIBRARIES 
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation 
    lardataobj::Simulation
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::AnalysisBase
    lardata::Utilities
    lardataalg::DetectorInfo
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    ROOT::Core
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    larreco::RecoAlg
    dunereco::AnaUtils
    dunereco::Ivysaurus_Utils
    larpandora::LArPandoraInterface
    larreco::Calorimetry
)

install_headers()
install_fhicl()
install_source()

