cmake_minimum_required(VERSION 3.16)
project(study1_app CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(Geant4 REQUIRED)  # conda geant4 exports Geant4::Geant4 aggregate
find_package(nlohmann_json QUIET)
if(NOT nlohmann_json_FOUND)
  # conda env ships nlohmann_json; fall back to header-only FetchContent
  include(FetchContent)
  FetchContent_Declare(json HEADER_ONLY
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
  FetchContent_MakeAvailable(json)
endif()
add_executable(study1_app src/StudyApp.cc)
# conda Geant4 exports per-library targets only (no Geant4::Geant4 aggregate)
target_link_libraries(study1_app PRIVATE
  Geant4::G4global Geant4::G4run Geant4::G4event Geant4::G4tracking
  Geant4::G4track Geant4::G4particles Geant4::G4geometry Geant4::G4materials
  Geant4::G4processes Geant4::G4processes_core Geant4::G4processes_hadronic
  Geant4::G4physicslists Geant4::G4intercoms Geant4::G4ptl
  Geant4::G4digits_hits Geant4::G4analysis Geant4::G4readout
  Geant4::G4vis_management Geant4::G4modeling Geant4::G4graphics_reps
  Geant4::G4interfaces Geant4::G4parmodels Geant4::G4error_propagation
  Geant4::G4mctruth Geant4::G4geomtext Geant4::G4tools
  nlohmann_json::nlohmann_json)
target_include_directories(study1_app PRIVATE src)