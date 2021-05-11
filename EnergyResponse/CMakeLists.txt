
# which cmake to use at least, otherwise fail
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# set up cmake project
project(E4SimulationAppProject)

# tell the compiler where to find the actual files
include_directories(${PROJECT_SOURCE_DIR}/include)

# collect names of the user code source files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

# tell the compiler which source files need to be compiled
add_executable(sim main.cc ${sources})

# locate Geant4 installation
# include packages for user interface (ui_all) and visualization (vis_all)
find_package(Geant4 REQUIRED ui_all vis_all)

# include compilation definitions for Geant4
include(${Geant4_USE_FILE})

# link against Geant4 libraries
target_link_libraries(sim ${Geant4_LIBRARIES})

# copy the init_vis macro to te build directory
set(SIM_SCRIPTS macros/init_vis.mac)

foreach(_script ${SIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
