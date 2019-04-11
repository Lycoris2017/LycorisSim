
#--------------------------------------------------------
# Setup CMake
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(LycorisSim)

SET(CMAKE_INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib")

#--------------------------------------------------------
# Setup ROOT, need for analysis code
#

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#--------------------------------------------------------
# Find G4
#

find_package(Geant4 REQUIRED)
message(STATUS "Geant4_CXXSTD=${Geant4_CXXSTD}")


#--------------------------------------------------------
# Setup G4
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#--------------------------------------------------------
# Setup project
#

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#--------------------------------------------------------
# Add Executables
#
add_executable(tpcsim main/tpcsim.cc ${sources} ${headers})
target_link_libraries(tpcsim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

add_executable(old_tpcsim main/old_tpcsim.cc ${sources} ${headers})
target_link_libraries(old_tpcsim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#--------------------------------------------------------
# Install Exe
#
install(TARGETS tpcsim old_tpcsim
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
