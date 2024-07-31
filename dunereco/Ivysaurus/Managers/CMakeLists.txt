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
    larreco::Calorimetry
    lardata::Utilities
    lardataalg::DetectorInfo
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    ROOT::Core
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    larpandora::LArPandoraInterface
    dunereco::AnaUtils
    larreco::RecoAlg
    dunereco::Ivysaurus_Utils
)

install_headers()
install_fhicl()
install_source()
install_scripts()
