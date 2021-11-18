project(demag_nonequi_layers_cpp)
cmake_minimum_required(VERSION 3.15)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # for IDE support

find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python numpy REQUIRED)

set(CMAKE_SHARED_MODULE_PREFIX "") # skip the 'lib' prefix: lib<x>.so -> <x>.so

set(LIB_DEMAG_NONEQUI_LAYERS ${PROJECT_NAME}) # use project name as libname

add_library(${LIB_DEMAG_NONEQUI_LAYERS} MODULE demag_nonequi_layers/lib/demag_nonequi_layers.cpp)

target_compile_features(${LIB_DEMAG_NONEQUI_LAYERS} PRIVATE cxx_std_11) # min required std
target_compile_options(${LIB_DEMAG_NONEQUI_LAYERS} PRIVATE -O3 -Wall -Wpedantic) # -Wconversion

target_link_libraries(${LIB_DEMAG_NONEQUI_LAYERS} ${Boost_LIBRARIES})
target_include_directories(${LIB_DEMAG_NONEQUI_LAYERS} PRIVATE ${PYTHON_INCLUDE_DIRS})

include(GNUInstallDirs)
install(TARGETS ${LIB_DEMAG_NONEQUI_LAYERS} LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})

option(PACK_AS_DEB "Pack as .deb file instead of package script and .tar files" ON)
if(PACK_AS_DEB)
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "pthon")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON) # resolve dependencies automatically
endif(PACK_AS_DEB)

include(CPack)
