cmake_minimum_required (VERSION 3.25)

option (BUILD_THIRD_PARTIES "Whether or not a superbuild should be invoked." OFF)
if(BUILD_THIRD_PARTIES)
    project(suprebuild NONE) 
    include (${CMAKE_CURRENT_LIST_DIR}/cmake/buildpackages.cmake)
    set(BUILD_THIRD_PARTIES OFF CACHE BOOL "Whether or not a superbuild should be invoked")
    return()
endif()

include (${CMAKE_CURRENT_LIST_DIR}/cmake/config.cmake)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

list(APPEND REQUIRED_FEATURES 
    cxx_std_20
)


project(AtlasSeg 
    VERSION 0.0.1
    LANGUAGES CXX
)


find_package(ITK REQUIRED)

if(NOT ITK_FOUND)
    message(FATAL_ERROR "ITK package not found!!")
else()
    message(STATUS "ITK package found at ${ITK_DIR}")
endif()


include(${ITK_USE_FILE})

set(CMAKE_SKIP_BUILD_RPATH TRUE)
# Note that shared libraries are install into /lib or /bin directory for linux or windows. 
# set(ITK_LIBRARY_PATH "${CMAKE_CURRENT_LIST_DIR}/thirdparties/itk/$<IF:$<BOOL:${WIN32}>,bin,lib>/")
set(ITK_LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/lib")

# set(ENV{LD_LIBRARY_PATH} "${ITK_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
# Specifies the runtime search path (rpath) for the executable. Multiple paths are supported.
set(CMAKE_INSTALL_RPATH ${ITK_LIBRARY_PATH})
# Build into the rpath into the executable, to avoid the library not found error when you run it.
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# Automatically adds directories that contain libraries linked to the executable to the rpath
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

set(PROJECT_SOVERSION 0)

SET(SOURCE_FILES 
    ${CMAKE_CURRENT_LIST_DIR}/src/spline.cpp
)

SET(INCLUDE_FILES
    # ${CMAKE_CURRENT_LIST_DIR}/inc
)

# Retrieve delta path of current dir from root dir
string(REPLACE ${CMAKE_SOURCE_DIR} "" TARGET_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(AtlasSeg  ${SOURCE_FILES})

set_target_properties(AtlasSeg PROPERTIES
    OUTPUT_NAME AtlasSeg
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_SOVERSION}
)
# set_target_properties(AtlasSeg PROPERTIES
#     INSTALL_RPATH ${ITK_LIBRARY_PATH}
#     BUILD_WITH_INSTALL_RPATH TRUE
# )



target_compile_features(AtlasSeg PUBLIC ${REQUIRED_FEATURES})

target_include_directories(AtlasSeg PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>                     # include folder 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/inc>
)

target_link_libraries(AtlasSeg PUBLIC ${ITK_LIBRARIES})


# Generates VS projects in a folder structure that mirrors that of the source code 
if(MSVC)
        
    if(NOT "${TARGET_FOLDER}" STREQUAL "")

        set_target_properties(${SHARED_LIB} PROPERTIES FOLDER ${TARGET_FOLDER})

    endif()

    source_group("src" FILES ${SOURCE_FILES})

    source_group("inc" FILES ${INCLUD_FILES})

    source_group("res" FILES ${RESOURCE_FILES})

endif()


install(TARGETS AtlasSeg
    # RUNTIME_DEPENDENCIES
    # DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/thirdparties/itk/lib"
    # PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
    # POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} #"$<IF:$<BOOL:${WIN32}>,bin,bin>" # dll
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # shared libraries but not dll (not useful)
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # mrobjects_static.lib
    OBJECTS DESTINATION ${CMAKE_INSTALL_LIBDIR} # installation for object files *.obj
    COMPONENT Runtime
)




enable_testing()

#add_subdirectory(tst)
#]==]