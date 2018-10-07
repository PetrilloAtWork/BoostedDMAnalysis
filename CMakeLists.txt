# ======================================================================
#  boosteddmanalysis main build file
# ======================================================================

cmake_minimum_required (VERSION 3.10)

project(boosteddmanalysis)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules
                      $ENV{CETBUILDTOOLS_DIR}/Modules
                      ${CMAKE_MODULE_PATH}
                      )

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS VIGILANT
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
  )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larsoft v1_00_00 )
find_ups_product( nutools v1_00_00 )
find_ups_product( art v2_10_00 )
find_ups_product( cetbuildtools v4_14_01 )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source code subdirectories
add_subdirectory(boosteddmanalysis)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
