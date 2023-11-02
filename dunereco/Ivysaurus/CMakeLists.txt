include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  add_subdirectory(Modules)
  add_subdirectory(Managers)
endif ()

install_headers()
install_fhicl()
install_source()

