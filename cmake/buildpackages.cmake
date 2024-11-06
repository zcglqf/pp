include(ExternalProject)

# find package

message(STATUS "Start building ITK ...")
ExternalProject_Add(
    ITK
    GIT_REPOSITORY https://github.com/InsightSoftwareConsortium/ITK.git
    GIT_TAG        v5.4.0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ITK
    CMAKE_ARGS 
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_STATIC_LIBS=ON
        -DBUILD_OBJECT_LIBS=ON
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_LIST_DIR}/../thirdparties/itk
    LOG_DOWNLOAD ON
    LOG_UPDATE ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
)
message(STATUS "End building ITK ...")

message(STATUS "Start building ANTS")
ExternalProject_Add(
  ANTs
  GIT_REPOSITORY https://github.com/ANTsX/ANTs.git
  GIT_TAG v2.5.3
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ANTs
  CMAKE_ARGS
      -DBUILD_SHARED_LIBS=ON 
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_LIST_DIR}/../thirdparties/ants
  DEPENDS ITK  # ANTs will depend on the ITK project
  
)
message(STATUS "End building ANTS ...")