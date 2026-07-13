cmake_minimum_required(VERSION 3.14)

project(
    Deexcitation
    VERSION 1.0.0
    LANGUAGES CXX
)

find_package(Geant4 REQUIRED)

find_package(COLA REQUIRED)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME Deexcitation)

file(GLOB SOURCES CONFIGURE_DEPENDS Deexcitation/*.cpp Deexcitation/handler/*.cpp)

add_library(${COLA_MODULE_NAME} SHARED ${SOURCES})

target_link_libraries(${COLA_MODULE_NAME} PUBLIC COLA ${Geant4_LIBRARIES})
target_include_directories(${COLA_MODULE_NAME} PUBLIC ${Geant4_INCLUDE_DIR})
target_compile_options(${COLA_MODULE_NAME} PRIVATE -Wall -Werror -Wextra -Wpedantic)
target_compile_features(${COLA_MODULE_NAME} PRIVATE cxx_std_20)

target_include_directories(
    ${COLA_MODULE_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Deexcitation>
)

set_target_properties(
    ${COLA_MODULE_NAME}
    PROPERTIES VERSION "${PROJECT_VERSION}"
               SOVERSION "${PROJECT_VERSION_MAJOR}"
               BUILD_WITH_INSTALL_RPATH TRUE
               INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Check for COLA_FETCH_MODE to prevent install instructions from generating
if (NOT DEFINED ENV{COLA_FETCH_MODE})
    install(
        TARGETS ${COLA_MODULE_NAME}
        EXPORT ${COLA_MODULE_NAME}Export
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${COLA_MODULE_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${COLA_MODULE_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
    )

    install(
        DIRECTORY ${COLA_MODULE_NAME}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
        PATTERN "*.h"
    )

    configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
    )

    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
        VERSION ${Deexcitation_VERSION}
        COMPATIBILITY AnyNewerVersion
    )

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
    )

    install(
        EXPORT ${COLA_MODULE_NAME}Export
        FILE ${COLA_MODULE_NAME}Export.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
    )

    enable_testing()
    if (BUILD_TESTING)
        add_subdirectory(tests)
    endif ()

    setup_cola_quality_tools(EXCLUDE_PATTERNS "(^|/)geant4/" "(^|/)tests/")
endif ()
