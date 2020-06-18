# The SUDODEM has the following parameters to configure:
#  INSTALL_PREFIX: path, where SudoDEM will be installed (/usr/local by default)
#  LIBRARY_OUTPUT_PATH: path to install libraries (lib by default)
#  DEBUG: compile in debug-mode (OFF by default)
#  CLUSTER: compile on a cluster
#  CMAKE_VERBOSE_MAKEFILE: output additional information during compiling (OFF by default)
#  SUFFIX: suffix, added after binary-names (version number by default)
#  NOSUFFIX: do not add a suffix after binary-name (OFF by default)
#  SUDODEM_VERSION: explicitely set version number (is defined from git-directory by default)
#  ENABLE_GUI: enable GUI option (ON by default)
#  ENABLE_OPENMP: enable OpenMP-parallelizing option (ON by default)
#  runtimePREFIX: used for packaging, when install directory is not the same is runtime directory (/usr/local by default)
#  CHUNKSIZE: set >1, if you want several sources to be compiled at once. Increases compilation speed and RAM-consumption during it (1 by default).
#  VECTORIZE: enables vectorization and alignment in Eigen3 library, experimental (OFF by default)

project(sudoDEM C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")
set(Python_ADDITIONAL_VERSIONS 2.7)###zhswee

find_program(LSB_RELEASE_EXEC lsb_release)
execute_process(COMMAND ${LSB_RELEASE_EXEC} -is
    OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

OPTION(CLUSTER "Compilation on a cluster" ON)
IF ( ${LSB_RELEASE_ID_SHORT} MATCHES "Ubuntu")
  SET(CLUSTER OFF)
ENDIF ( ${LSB_RELEASE_ID_SHORT} MATCHES "Ubuntu")

MESSAGE("The distribution is: ${LSB_RELEASE_ID_SHORT}; and Option CLUSTER is set as ${CLUSTER}")

#INCLUDE(FindPythonLibs)
#INCLUDE(FindPythonInterp)
find_package( PythonInterp 2.7 REQUIRED )

IF (CLUSTER)
#find_package( PythonLibs 2.7 REQUIRED PATHS '/opt/rh/python27/root/usr/lib64')
MESSAGE("Specify Python lib on a cluster manually.")
set(PythonLibs_DIR '/opt/rh/python27/root/usr/lib64')
set(PYTHON_LIBRARIES /opt/rh/python27/root/usr/lib64/libpython2.7.so)
ELSE (CLUSTER)
find_package( PythonLibs 2.7 REQUIRED)
MESSAGE("Specify Python lib on a cluster automatically.")
ENDIF (CLUSTER)

#find_package( PythonLibs )

INCLUDE(FindOpenMP)
INCLUDE(FindQt4)
INCLUDE(FindPkgConfig)
INCLUDE(GetVersion)
INCLUDE(FindOpenGL)
INCLUDE(FindNumPy)

INCLUDE(FindPythonModule)
INCLUDE(GNUInstallDirs)

MESSAGE(STATUS "PYTHON_LIBRARIES:${PYTHON_LIBRARIES}")

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
#===========================================================
# Add possibility to use local boost installation (e.g. -DLocalBoost=1.46.1)
#set(CMAKE_MODULE_PATH "/scratch/softwares/boostinstall/lib/;${CMAKE_MODULE_PATH}")
set(Boost_NO_SYSTEM_PATHS TRUE)
#if (Boost_NO_SYSTEM_PATHS)
IF (CLUSTER)
  set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../libsinstall")
ELSE (CLUSTER)
  set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib")
ENDIF (CLUSTER)
  #set(BOOST_ROOT "/tmp/project/")
  set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
  set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
  #set(BOOST_LIBRARYDIR "${BOOST_ROOT}/lib")
  #set(BOOST_LIBRARYDIR "/opt/SudoDEM/lib/3rdlibs/")


#endif (Boost_NO_SYSTEM_PATHS)
FIND_PACKAGE(Boost 1.67 COMPONENTS python27 thread filesystem iostreams regex serialization system date_time REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS python thread filesystem iostreams regex serialization system date_time REQUIRED PATHS ${BOOST_LIBRARYDIR} NO_DEFAULT_PATH)

#set(Boost_LIBRARIES "-lboost_python -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_regex -lboost_serialization -lboost_system -lboost_date_time -lpthread") #"/home/sway/software/DEM/3rdlib/boost159/lib/libboost_python.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_thread.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_filesystem.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_iostreams.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_regex.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_serialization.so/home/sway/software/DEM/3rdlib/boost159/lib/libboost_system.so/home/sway/software/DEM/3rdlib/boost159/lib/libboos#t_date_time.so/usr/lib/x86_64-linux-gnu/libpthread.so")
#set(Boost_LIBRARIES "${BOOST_LIBRARY_DIRS}libboost_python.so${BOOST_LIBRARY_DIRS}libboost_thread.so${BOOST_LIBRARY_DIRS}libboost_filesystem.so${BOOST_LIBRARY_DIRS}libboost_iostreams.so${BOOST_LIBRARY_DIRS}libboost_regex.so${BOOST_LIBRARY_DIRS}libboost_serialization.so${BOOST_LIBRARY_DIRS}libboost_system.so${BOOST_LIBRARY_DIRS}libbost_date_time.so")
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
#LINK_DIRECTORIES( ${BOOST_LIBRARYDIR} )
# for checking purpose
MESSAGE("--   Boost_VERSION: " ${Boost_VERSION})
MESSAGE("--   Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
MESSAGE("--   BOOST_LIBRARYDIR: " ${BOOST_LIBRARYDIR})
MESSAGE("--   Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
MESSAGE("--   Boost_LIBRARIES: " ${Boost_LIBRARIES})

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
OPTION(ENABLE_VTK "Enable VTK" ${DEFAULT})
OPTION(ENABLE_OPENMP "Enable OpenMP" ${DEFAULT})
OPTION(ENABLE_GUI "Enable GUI" ${DEFAULT})

IF (CLUSTER)
SET (ENABLE_GUI OFF)
ENDIF (CLUSTER)
#===========================================================
# Use Eigen3 by default
# Use Eigen3 by default
#IF (EIGEN3_FOUND)
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/")
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  MESSAGE(STATUS "Found Eigen3, at: ${EIGEN3_INCLUDE_DIR}")

  # Minimal required version 3.2.1
  #IF ((${EIGEN3_MAJOR_VERSION} LESS 2) OR ((${EIGEN3_MAJOR_VERSION} EQUAL 2) AND (${EIGEN3_MINOR_VERSION} LESS 1)))
 #   MESSAGE(FATAL_ERROR "Minimal required Eigen3 version is 3.2.1, please update Eigen3!")
 # ENDIF ((${EIGEN3_MAJOR_VERSION} LESS 2) OR ((${EIGEN3_MAJOR_VERSION} EQUAL 2) AND (${EIGEN3_MINOR_VERSION} LESS 1)))

  IF (NOT VECTORIZE)
    MESSAGE(STATUS "Disable vectorization")
    ADD_DEFINITIONS("-DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT")
  ELSE (NOT VECTORIZE)
    MESSAGE(STATUS "Enable vectorization")
  ENDIF (NOT VECTORIZE)

#ENDIF(EIGEN3_FOUND)
#===========================================================
INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
SET(LINKLIBS  "${LINKLIBS};${BZIP2_LIBRARIES};${ZLIB_LIBRARIES};")

SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} Odeint")
ADD_DEFINITIONS("-DSUDODEM_ODEINT")
#===========================================================
IF(ENABLE_OPENMP)
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
ELSE(ENABLE_OPENMP)
  SET(DISABLED_FEATS "${DISABLED_FEATS} OPENMP")
ENDIF(ENABLE_OPENMP)
#===========================================================
IF(ENABLE_GUI)
  FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL)
  FIND_PACKAGE(OpenGL)
  FIND_PACKAGE(GLUT)

  FIND_PACKAGE(glib2)
  #FIND_PACKAGE(QGLVIEWER)
  set(QGLVIEWER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib")
 set(QGLVIEWER_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/HeaderLib/lib/libQGLViewer.so")
 #set(QGLVIEWER_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/../3rdlib/libQGLViewer-2.6.3/QGLViewer/libQGLViewer.a")#using static library
  SET(QGLVIEWER_FOUND TRUE)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQGLVIEWER_FOUND")
SET(LINKLIBS  "${LINKLIBS};${QGLVIEWER_LIBRARIES}")

  IF(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    SET(GUI_LIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} ${QGLVIEWER_LIBRARIES})
    SET(GUI_SRC_LIB "lib/opengl/GLUtils.cpp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSUDODEM_OPENGL")
    INCLUDE_DIRECTORIES(${GLIB2_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${QT_INCLUDES})

    MESSAGE(STATUS "Found GUI-LIBS")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} GUI")
  ELSE(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    MESSAGE(STATUS "GUI-LIBS NOT found")
    MESSAGE(STATUS "qt4 ${QT4_FOUND} opengl${OPENGL_FOUND} glut  ${GLUT_FOUND} glib2 ${GLIB2_FOUND} qglviewer ${QGLVIEWER_FOUND}")
    SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
    SET(ENABLE_GUI OFF)
  ENDIF(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
ELSE(ENABLE_GUI)
  SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
ENDIF(ENABLE_GUI)

#===========================================================

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#===========================================================

IF (NOT INSTALL_PREFIX)
 #SET(CMAKE_INSTALL_PREFIX "/usr/local")
 IF (CLUSTER)
 SET(CMAKE_INSTALL_PREFIX "~/software/FDEMpackage/SudoDEM")
 ELSE (CLUSTER)
 SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../sudodeminstall/SudoDEM2D")
 ENDIF (CLUSTER)
 MESSAGE("SudoDEM will be installed to default path ${CMAKE_INSTALL_PREFIX}, if you want to override it use -DINSTALL_PREFIX option.")
ELSE (NOT INSTALL_PREFIX)
	GET_FILENAME_COMPONENT(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} ABSOLUTE)
	MESSAGE("SudoDEM will be installed to ${CMAKE_INSTALL_PREFIX}")
ENDIF (NOT INSTALL_PREFIX)

IF (NOT SUFFIX)
  SET (SUFFIX "-${SUDODEM_VERSION}")
  #SET (SUFFIX "")
ENDIF (NOT SUFFIX)

IF(NOSUFFIX)   #For packaging
  SET (SUFFIX "")
ENDIF(NOSUFFIX)   #For packaging
SET (LIBRARY_OUTPUT_PATH lib)
IF(NOT LIBRARY_OUTPUT_PATH)   #For packaging
  SET (LIBRARY_OUTPUT_PATH ${CMAKE_INSTALL_LIBDIR})
ENDIF(NOT LIBRARY_OUTPUT_PATH)   #For packaging

IF (NOT runtimePREFIX)
  SET (runtimePREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF (NOT runtimePREFIX)

MESSAGE (STATUS "Suffix is set to " ${SUFFIX})
MESSAGE (STATUS "LIBRARY_OUTPUT_PATH is set to " ${LIBRARY_OUTPUT_PATH})
MESSAGE (STATUS "runtimePREFIX is set to " ${runtimePREFIX})
#===========================================================

#SET(SUDODEM_LIB_PATH ${CMAKE_INSTALL_PREFIX}/${LIBRARY_OUTPUT_PATH}/sudodem${SUFFIX})
SET(SUDODEM_LIB_PATH ${CMAKE_INSTALL_PREFIX}/lib/sudodem${SUFFIX})
SET(SUDODEM_EXEC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
SET(SUDODEM_PY_PATH ${SUDODEM_LIB_PATH}/py)
SET(SUDODEM_DOC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/doc/sudodem${SUFFIX})
SET(SUDODEM_MAN_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_MANDIR})


SET(CMAKE_SKIP_BUILD_RPATH FALSE)

IF (CLUSTER)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ELSE (CLUSTER)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  #SET(CMAKE_INSTALL_RPATH "${SUDODEM_LIB_PATH};${SUDODEM_PY_PATH};${SUDODEM_PY_PATH}/sudodem/;${SUDODEM_PY_PATH}/sudodem/qt;${QGLVIEWER_INCLUDE_DIR}/QGLViewer")
  #SET(CMAKE_INSTALL_RPATH "${SUDODEM_LIB_PATH};${SUDODEM_PY_PATH};${SUDODEM_PY_PATH}/sudodem/;${SUDODEM_PY_PATH}/sudodem/qt;$ORIGIN/../../../../../3rdlibs")
  SET(CMAKE_INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")#this is specified to python modules in ${SUDODEM_PY_PATH}/sudodem/.

  #SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
ENDIF (CLUSTER)

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

#SET(QHULL_INCLUDE_DIR,"lib/")
#MESSAGE(STATUS "Found QHULL: ${QHULL_INCLUDE_DIR}, ${QHULL_LIBRARIES}")
#INCLUDE_DIRECTORIES(${QHULL_INCLUDE_DIR})

FILE(GLOB QHULL_SRC_LIB "lib/qhull/*.c")
ADD_LIBRARY(qhull SHARED ${QHULL_SRC_LIB})
SET_TARGET_PROPERTIES(qhull PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../3rdlibs/")
TARGET_LINK_LIBRARIES(qhull)
INSTALL(TARGETS qhull DESTINATION ${SUDODEM_LIB_PATH})
#===========================================================
FILE(GLOB VORO_SRC_LIB "lib/voro++/voro++.cc")
ADD_LIBRARY(voro++ SHARED ${VORO_SRC_LIB})
SET_TARGET_PROPERTIES(voro++ PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../3rdlibs/")
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
IF (NOT CLUSTER)
find_python_module(Tkinter REQUIRED)
ENDIF (NOT CLUSTER)
#===========================================================


ADD_EXECUTABLE(sudodem2d ${CMAKE_CURRENT_SOURCE_DIR}/core/main/sudodem.cpp)
SET_TARGET_PROPERTIES(sudodem2d PROPERTIES INSTALL_RPATH "$ORIGIN;$ORIGIN/../lib/3rdlibs/;")

ADD_LIBRARY(boot SHARED ${CMAKE_CURRENT_SOURCE_DIR}/core/main/pyboot.cpp)
SET_TARGET_PROPERTIES(boot PROPERTIES PREFIX "" LINK_FLAGS "-Wl,-z,origin,--as-needed" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs;") #add -z,origin
TARGET_LINK_LIBRARIES(sudodem qhull ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${LINKLIBS} -lrt)

MESSAGE("link libraries: " ${LINKLIBS})
#TARGET_LINK_LIBRARIES(sudodem qhull -lboost_python -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_regex -lboost_serialization -lboost_system -lboost_date_time -lpthread ${PYTHON_LIBRARIES} ${LINKLIBS} -lrt)
#TARGET_LINK_LIBRARIES(sudodem qhull -L/opt/SudoDEM/lib/3rdlibslibboost_python.so libboost_thread.so libboost_filesystem.so libboost_iostreams.so libboost_regex.so libboost_serialization.so libboost_system.so libboost_date_time.so libthread.so ${PYTHON_LIBRARIES} ${LINKLIBS} -lrt)
#SET_TARGET_PROPERTIES(sudodem  PROPERTIES LINK_FLAGS "-Wl,-z,origin,--as-needed" ) #add -z,origin
SET_TARGET_PROPERTIES(sudodem  PROPERTIES LINK_FLAGS "-Wl,--as-needed" INSTALL_RPATH "$ORIGIN;$ORIGIN/../3rdlibs;$ORIGIN/../3rdlibs/py;$ORIGIN/py/sudodem/qt" )
TARGET_LINK_LIBRARIES(boot sudodem)

IF (CLUSTER)
  TARGET_LINK_LIBRARIES(sudodem2d ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
ELSE (CLUSTER)
  TARGET_LINK_LIBRARIES(sudodem2d boot ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})#FIXME:sudodem2d is not dependent on boot on Ubuntu, but not on CentOS
ENDIF (CLUSTER)

IF(ENABLE_GUI)
  TARGET_LINK_LIBRARIES(sudodem _GLViewer ${GUI_LIBS})
ENDIF(ENABLE_GUI)

#====================================
#IF (NOT (PY_minieigen))
#  ADD_LIBRARY(miniEigen SHARED py/mathWrap/miniEigen.cpp)
#  SET_TARGET_PROPERTIES(miniEigen PROPERTIES PREFIX "")
#  TARGET_LINK_LIBRARIES(miniEigen ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
# INSTALL(TARGETS miniEigen DESTINATION ${SUDODEM_PY_PATH})
#ENDIF (NOT (PY_minieigen))

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
