cmake_minimum_required(VERSION 3.16)

project(SimplexSolver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets PrintSupport)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS
    Widgets
    PrintSupport
)

set(PROJECT_SOURCES
    main.cpp
    mainwindow.cpp
    mainwindow.ui
    modelo.cpp
    simplexsolver.cpp
    uibuilder.cpp
    qcustomplot.cpp
    resources.qrc
    app_icon.rc
)

set(PROJECT_HEADERS
    mainwindow.h
    modelo.h
    simplexsolver.h
    uibuilder.h
    qcustomplot.h
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(SimplexSolver
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
        ${PROJECT_HEADERS}
        resources.qrc
    )
else()
    if(ANDROID)
        add_library(SimplexSolver SHARED
            ${PROJECT_SOURCES}
            ${PROJECT_HEADERS}
        )
    else()
        add_executable(SimplexSolver
            ${PROJECT_SOURCES}
            ${PROJECT_HEADERS}
        )
    endif()
endif()

target_link_libraries(SimplexSolver PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::PrintSupport
)

if(${QT_VERSION} VERSION_LESS 6.1.0)
    set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.yourorganization.SimplexSolver)
endif()

set_target_properties(SimplexSolver PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS SimplexSolver
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION bin
)

install(FILES "${CMAKE_SOURCE_DIR}/styles.qss" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/grafico3D.py" DESTINATION .)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/assets/" DESTINATION .)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/_python_env/" DESTINATION .)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(SimplexSolver)
endif()

include(CPack)

set(CPACK_PACKAGE_NAME "SimplexSolver")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_VENDOR "Claudia Rocha")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A desktop application to solve Linear Programming problems with 3D visualization.")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/Claudia-Rocha-H/SimplexSolverQt/")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")

if(WIN32)
    set(CPACK_GENERATOR "WIX")
    set(CPACK_WIX_EXTENSIONS "WixUtilExtension;WixUIExtension")
    set(CPACK_WIX_DEPLOY_QT ON)
    set(CPACK_WIX_SKIP_INSTALL_DIRECTORIES ON)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-win64")

elseif(APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_MAC_BUNDLE_CONTENT_DIRS "bin;Resources")
    set(CPACK_MAC_DEPLOY_QT ON)
    set(CPACK_BUNDLE_GUI_NAME "Simplex Solver")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-macos")

elseif(UNIX)
    set(CPACK_GENERATOR "DEB;RPM;TGZ")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libqt6core6, libqt6gui6, libqt6widgets6")
    set(CPACK_RPM_PACKAGE_REQUIRES "qt6-qtbase")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-linux")
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
