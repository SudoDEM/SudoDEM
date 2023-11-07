project(sudoDEM C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")
set(Python_ADDITIONAL_VERSIONS 3.8)
find_package( PythonInterp 3.8 REQUIRED )
find_package( PythonLibs 3.8 REQUIRED)

INCLUDE(FindOpenMP)
INCLUDE(FindPkgConfig)
INCLUDE(GetVersion)
INCLUDE(FindOpenGL)
INCLUDE(FindNumPy)

INCLUDE(FindPythonModule)
INCLUDE(GNUInstallDirs)


#===========================================================
# HACK!!! If the version of gcc is 4.8 or greater, we add -ftrack-macro-expansion=0
# and -save-temps into compiler to reduce the memory consumption during compilation.
# See http://bugs.debian.org/726009 for more information
# Can be removed later, if gcc fixes its regression
# Taken from http://stackoverflow.com/questions/4058565/check-gcc-minor-in-cmake

EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
IF (GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
  MESSAGE(STATUS "GCC Version >= 4.8. Adding -ftrack-macro-expansion=0 and -save-temps")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0 -save-temps")
ENDIF()

#===========================================================

IF ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*clang")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -ftemplate-depth-512")
ENDIF()

#===========================================================

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")#sway

#===========================================================

IF (CMAKE_CXX_FLAGS)
  #If flags are set, add only neccessary flags
  IF (DEBUG)
    SET(CMAKE_VERBOSE_MAKEFILE 1)
    SET(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS("-DSUDODEM_DEBUG")
  ELSE (DEBUG)
    SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")
  ENDIF (DEBUG)
ELSE (CMAKE_CXX_FLAGS)
  #If flags are not set, add all useful flags
  IF (DEBUG)
    SET(CMAKE_VERBOSE_MAKEFILE 1)
    SET(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS("-DSUDODEM_DEBUG")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fPIC -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
  ELSE (DEBUG)
    SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS  "-Wall -fPIC -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -s")
  ENDIF (DEBUG)
ENDIF (CMAKE_CXX_FLAGS)

#=======================  BOOST LIBRARY  =======================#  
set(BOOST_NO_SYSTEM_PATHS TRUE)
set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib")
set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
set(BOOST_LIBRARIES "-lboost_python38 -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_regex -lboost_serialization -lboost_system -lboost_date_time -lpthread")
INCLUDE_DIRECTORIES (${BOOST_INCLUDE_DIRS})
LINK_DIRECTORIES( ${BOOST_LIBRARY_DIRS} )
MESSAGE(STATUS "BOOST_LIBRARYDIR: " ${BOOST_LIBRARY_DIRS})
MESSAGE(STATUS "BOOST_INCLUDE_DIRS: " ${BOOST_INCLUDE_DIRS})
MESSAGE(STATUS "BOOST_LIBRARIES: " ${BOOST_LIBRARIES})

#===========================================================
FIND_PACKAGE(NumPy REQUIRED)
INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIRS})

FIND_PACKAGE(BZip2 REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)
#===========================================================

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
SET(LINKLIBS "")
SET(CONFIGURED_FEATS "")
SET(DISABLED_FEATS "")

OPTION(NOSUFFIX "NOSUFFIX" ON)
OPTION(ENABLE_OPENMP "Enable OpenMP" ${DEFAULT})
OPTION(ENABLE_GUI "Enable GUI" ${DEFAULT})

#=======================  EIGEN  =======================# 
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/Eigen-3.3.5/")
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
MESSAGE(STATUS "Found Eigen3, at: ${EIGEN3_INCLUDE_DIR}")

IF (NOT VECTORIZE)
  MESSAGE(STATUS "Disable vectorization")
  ADD_DEFINITIONS("-DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT")
ELSE (NOT VECTORIZE)
  MESSAGE(STATUS "Enable vectorization")
ENDIF (NOT VECTORIZE)

#===========================================================
INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
SET(LINKLIBS  "${LINKLIBS};${BZIP2_LIBRARIES};${ZLIB_LIBRARIES};")

#===========================================================
SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} Odeint")
ADD_DEFINITIONS("-DSUDODEM_ODEINT")

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSUDODEM_OPENMP ${OpenMP_CXX_FLAGS}")
    MESSAGE(STATUS "Found OpenMP")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} OpenMP")
ELSE(OPENMP_FOUND)
    MESSAGE(STATUS "OpenMP NOT found")
    SET(ENABLE_OPENMP OFF)
    SET(DISABLED_FEATS "${DISABLED_FEATS} OPENMP")
ENDIF(OPENMP_FOUND)

#===========================================================

#=======================  QT5  ====================================
MESSAGE(STATUS "Enable GUI: " ${ENABLE_GUI})
IF(ENABLE_GUI)
  FIND_PACKAGE(Qt5 CONFIG REQUIRED Widgets Xml OpenGL)
  SET(QGLVIEWER_FOUND TRUE)
  FIND_PACKAGE(Qt5Widgets)
  FIND_PACKAGE(OpenGL)
  FIND_PACKAGE(GLUT)
  FIND_PACKAGE(glib2)


  set(QGLVIEWER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib")
  INCLUDE_DIRECTORIES(${QGLVIEWER_INCLUDE_DIR})
  set(QGLVIEWER_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib/lib/libQGLViewer.so")
  SET(QGLVIEWER_FOUND TRUE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQGLVIEWER_FOUND")
  SET(LINKLIBS  "${LINKLIBS};${QGLVIEWER_LIBRARIES}")

  IF(Qt5Widgets_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    SET(GUI_LIBS ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${QGLVIEWER_LIBRARIES})
    SET(GUI_SRC_LIB "lib/opengl/GLUtils.cpp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSUDODEM_OPENGL")
    #INCLUDE_DIRECTORIES(${GLIB2_INCLUDE_DIRS})
    #INCLUDE_DIRECTORIES(${QT_INCLUDES})

    MESSAGE(STATUS "Found GUI-LIBS")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} GUI")
  ELSE(Qt5Widgets_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    MESSAGE(STATUS "GUI-LIBS NOT found")
    MESSAGE(STATUS "qt5 ${QT5Widgets_FOUND} opengl ${OPENGL_FOUND} glut  ${GLUT_FOUND} glib2 ${GLIB2_FOUND} qglviewer ${QGLVIEWER_FOUND}")
    SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
    SET(ENABLE_GUI OFF)
  ENDIF(Qt5Widgets_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
ELSE(ENABLE_GUI)
  SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
ENDIF(ENABLE_GUI)
#===============================================

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#===========================================================
# floating_point_utilities_v3 are already in Boost included
# Use embedded copy only if Boost older than 1.47.0

#IF((Boost_MINOR_VERSION LESS 47) AND (Boost_MAJOR_VERSION EQUAL 1))
#  MESSAGE(STATUS "Boost version is less than 1.47, using embedded version of #floating_point_utilities_v3")
#  MESSAGE(STATUS "Consider updating boost or system, as this embedded library will be removed soon")
#  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/extra/floating_point_utilities_v3)
#ENDIF((Boost_MINOR_VERSION LESS 47) AND (Boost_MAJOR_VERSION EQUAL 1))

#===========================================================
SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/../sudodeminstall/SudoDEM2D")
MESSAGE(STATUS "SudoDEM will be installed to ${CMAKE_INSTALL_PREFIX}")

#===========================================================
SET(SUDODEM_LIB_PATH ${CMAKE_INSTALL_PREFIX}/lib/sudodem)
SET(SUDODEM_EXEC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
SET(SUDODEM_PY_PATH ${SUDODEM_LIB_PATH}/py)
SET(SUDODEM_DOC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/doc/sudodem${SUFFIX})
SET(SUDODEM_MAN_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_MANDIR})


SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


#===========================================================
IF(ENABLE_GUI)
    ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/gui")
ENDIF(ENABLE_GUI)
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/py")
#===========================================================

# this is only a temporary hack, headers should be installed in the build-dir directly
# perhaps there is a proper way, have to ask
IF(NOT EXISTS "${CMAKE_BINARY_DIR}/sudodem")
  EXECUTE_PROCESS(COMMAND ln -s ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}/sudodem)
ENDIF ()
#===========================================================


FILE(GLOB QHULL_SRC_LIB "lib/qhull/*.c")
ADD_LIBRARY(qhull SHARED ${QHULL_SRC_LIB})
SET_TARGET_PROPERTIES(qhull PROPERTIES PREFIX "")
TARGET_LINK_LIBRARIES(qhull)
INSTALL(TARGETS qhull DESTINATION ${SUDODEM_LIB_PATH})
#===========================================================
FILE(GLOB VORO_SRC_LIB "lib/voro++/voro++.cc")
ADD_LIBRARY(voro++ SHARED ${VORO_SRC_LIB})
SET_TARGET_PROPERTIES(voro++ PROPERTIES PREFIX "")
TARGET_LINK_LIBRARIES(voro++)
INSTALL(TARGETS voro++ DESTINATION ${SUDODEM_LIB_PATH})
#===========================================================


FILE(GLOB SRC_CORE "core/*.cpp")
FILE(GLOB_RECURSE SRC_PKG  "pkg/*.cpp")
FILE(GLOB SRC_LIB  "lib/*.cpp")

SET(SRC_LIB "${SRC_LIB};lib/base/Math.cpp;lib/factory/ClassFactory.cpp;lib/factory/DynLibManager.cpp")
SET(SRC_LIB "${SRC_LIB};lib/serialization/Serializable.cpp;lib/pyutil/gil.cpp;core/main/pyboot.cpp;${GUI_SRC_LIB};${CGAL_SRC_LIB}")

#===========================================================

IF (CHUNKSIZE)
  INCLUDE(CombineSources)
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/core "${SRC_CORE}" ${CHUNKSIZE})
  FILE(GLOB SRC_CORE_COMBINED "${CMAKE_BINARY_DIR}/core.*.cpp")
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/pkg "${SRC_PKG}" ${CHUNKSIZE})
  FILE(GLOB SRC_PKG_COMBINED "${CMAKE_BINARY_DIR}/pkg.*.cpp")
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/lib "${SRC_LIB}" ${CHUNKSIZE})
  FILE(GLOB SRC_LIB_COMBINED "${CMAKE_BINARY_DIR}/lib.*.cpp")
  ADD_LIBRARY(sudodem SHARED ${SRC_LIB_COMBINED} ${SRC_CORE_COMBINED} ${SRC_PKG_COMBINED})
ELSE (CHUNKSIZE)
  ADD_LIBRARY(sudodem SHARED ${SRC_CORE} ${SRC_PKG} ${SRC_LIB})
ENDIF (CHUNKSIZE)

#===========================================================
find_python_module(minieigen)
IF (PY_minieigen)
  MESSAGE(STATUS "Use system minieigen version")
ELSE (PY_minieigen)
  MESSAGE(STATUS "Use embedded version of minieigen. Please, consider installing the corresponding package")
ENDIF (PY_minieigen)

find_python_module(tkinter REQUIRED)


#===========================================================


ADD_EXECUTABLE(sudodem2d ${CMAKE_CURRENT_SOURCE_DIR}/core/main/sudodem.cpp)
SET_TARGET_PROPERTIES(sudodem2d PROPERTIES INSTALL_RPATH "$ORIGIN;$ORIGIN/../lib/sudodem/py/sudodem;$ORIGIN/../lib/3rdlibs;")

ADD_LIBRARY(boot SHARED ${CMAKE_CURRENT_SOURCE_DIR}/core/main/pyboot.cpp)
SET_TARGET_PROPERTIES(boot PROPERTIES PREFIX "" LINK_FLAGS "-Wl,-z,origin,--as-needed" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs;") #add -z,origin

TARGET_LINK_LIBRARIES(sudodem qhull ${BOOST_LIBRARIES} ${PYTHON_LIBRARIES} ${LINKLIBS} -lrt)


SET_TARGET_PROPERTIES(sudodem  PROPERTIES LINK_FLAGS "-Wl,--as-needed" INSTALL_RPATH "$ORIGIN;$ORIGIN/../3rdlibs;$ORIGIN/../3rdlibs/py;$ORIGIN/py/sudodem/qt" )
TARGET_LINK_LIBRARIES(boot sudodem)

TARGET_LINK_LIBRARIES(sudodem2d boot ${BOOST_LIBRARIES} ${PYTHON_LIBRARIES})

IF(ENABLE_GUI)
  TARGET_LINK_LIBRARIES(sudodem _GLViewer ${GUI_LIBS})
  
ENDIF(ENABLE_GUI)


#====================================
#Back compatibility with scons
SET (realVersion ${SUDODEM_VERSION})
SET (version ${SUDODEM_VERSION})
SET (pyExecutable ${PYTHON_EXECUTABLE})
SET (profile "default")
SET (sourceRoot "${CMAKE_CURRENT_SOURCE_DIR}")
#====================================
CONFIGURE_FILE(core/main/sudodemcfg.h.in "${CMAKE_BINARY_DIR}/sudodemcfg.h")
CONFIGURE_FILE(py/__init__.py.in "${CMAKE_BINARY_DIR}/__init__.py")
#===========================================================

INSTALL(FILES "${CMAKE_BINARY_DIR}/__init__.py" DESTINATION ${SUDODEM_PY_PATH}/sudodem/)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/doc/sudodem-logo-note.png" DESTINATION "${SUDODEM_DOC_PATH}/img")

INSTALL(TARGETS boot DESTINATION "${SUDODEM_PY_PATH}/sudodem/")
INSTALL(TARGETS sudodem DESTINATION ${SUDODEM_LIB_PATH})
INSTALL(TARGETS sudodem2d DESTINATION ${SUDODEM_EXEC_PATH}/)

#===========================================================
MESSAGE(STATUS "===========================================================")
MESSAGE(STATUS "SudoDEM configured with following features:${CONFIGURED_FEATS}")
MESSAGE(STATUS "Disabled features:${DISABLED_FEATS}")
IF (DEBUG)
  MESSAGE(STATUS "Debug build")
  SET (debugbuild " (debug build)")
ELSE (DEBUG)
  MESSAGE(STATUS "Optimized build")
ENDIF (DEBUG)
IF (CHUNKSIZE)
  MESSAGE(STATUS "CHUNKSIZE is set to " ${CHUNKSIZE})
ENDIF (CHUNKSIZE)
MESSAGE(STATUS "===========================================================")
#===========================================================
