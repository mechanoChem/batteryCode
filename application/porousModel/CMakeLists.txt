##
#  CMake script
##

# Set the name of the project and target:
SET(TARGET "main")
#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(ENV{DEAL_II_DIR} "/nfs/mcfs_comp/home/rudraa/software/dealii-8.4.1/installWithTrilinos")
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#DEAL_II_INITIALIZE_CACHED_VARIABLES()
#PROJECT(${TARGET})
#DEAL_II_INVOKE_AUTOPILOT()



SET(src_DIR "/Users/wzhenlin/workspace/GitLab/batteryCode/src")
SET(inc_DIR "/Users/wzhenlin/workspace/GitLab/batteryCode/include")

FILE (GLOB TARGET_SRC "${src_DIR}/battery_base/*.cc" "${src_DIR}/battery_base/porousModel/*.cc" "${src_DIR}/ElectricChemoExpression/*.cc" "${src_DIR}/ElectricChemoExpression/porousModel/*.cc" "${src_DIR}/model/*.cc" "${src_DIR}/model/porousModel/*.cc" "${src_DIR}/supplementary/*.cc"  "src/*.cc" "electricChemoFormula/*.cc" "main.cc")


DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

INCLUDE_DIRECTORIES(${inc_DIR})

