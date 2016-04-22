cmake_minimum_required (VERSION 3.3)
project(spectrum C CXX Fortran)

set(CMAKE_VERBOSE_MAKEFILE off)

#add_compile_options(-fno-align-commons)

add_library(glow egrid.f ephoto.f maxt.f rcolum.f etrans.f exsect.f fieldm.f vquart.f gchem.f geomag.f solzen.f glow.f iri90.f nrlmsise00.f qback.f rout.f snoem.f snoemint.f ssflux.f)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)
#------ test program: Auroral emissions ----------
add_executable(aurexample aurexample.f)
target_link_libraries(aurexample glow)
#------ test program: high energy precipitation -------------
add_executable(hexexample hexexample.f)
target_link_libraries(hexexample glow)

