include(ExternalProject)

# find package
if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/../thirdparties/itk/bin/none")
    message(STATUS "Building itk ...")
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
else()
    message(STATUS "itk already built ...")
endif()

