#==================miniEigen=====================================
#find_python_module(minieigen)
#IF (PY_minieigen)
#  MESSAGE(STATUS "Use system minieigen version")
#ELSE (PY_minieigen)
#  MESSAGE(STATUS "Use embedded version of minieigen. Please, consider installing the corresponding package")
#  FILE(GLOB_RECURSE SRC_minieigen "${CMAKE_CURRENT_SOURCE_DIR}/mathWrap/*.cc")
#  ADD_LIBRARY(minieigen SHARED ${SRC_minieigen})
#  SET_TARGET_PROPERTIES(minieigen PROPERTIES PREFIX "")
#  TARGET_LINK_LIBRARIES(minieigen ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
#  INSTALL(TARGETS minieigen DESTINATION ${SUDODEM_PY_PATH})
#ENDIF (PY_minieigen)

#==================miniEigen=====================================

#==================pygts=========================================
IF(ENABLE_GTS)
  find_python_module(gts)
  IF (PY_gts)
    MESSAGE(STATUS "Use system gts version")
  ELSE (PY_gts)
    MESSAGE(STATUS "Use embedded version of gts. Please, consider installing the corresponding package")
    FILE(GLOB SRC_PYGTS "${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/pygts-0.3.1/*.c")
    ADD_LIBRARY(_gts SHARED ${SRC_PYGTS})
    SET_TARGET_PROPERTIES(_gts PROPERTIES PREFIX "")
    TARGET_LINK_LIBRARIES(_gts gts)
    INSTALL(FILES 3rd-party/pygts-0.3.1/__init__.py DESTINATION ${SUDODEM_LIB_PATH}/py/gts)
    INSTALL(FILES 3rd-party/pygts-0.3.1/pygts.py DESTINATION ${SUDODEM_LIB_PATH}/py/gts)
    INSTALL(TARGETS _gts DESTINATION ${SUDODEM_LIB_PATH}/py/gts)
  ENDIF (PY_gts)
ENDIF(ENABLE_GTS)
#==================pygts=========================================

#==================pyModules=========================================

FILE(GLOB filesPY "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
#SET(filesPY "${filesPY};${CMAKE_CURRENT_SOURCE_DIR}/pack/pack.py")
INSTALL(FILES ${filesPY} DESTINATION ${SUDODEM_PY_PATH}/sudodem)
INSTALL(FILES 3rd-party/mtTkinter-0.4/mtTkinter.py DESTINATION ${SUDODEM_LIB_PATH}/py)
#FILE(GLOB filesPYTests "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.py")
#INSTALL(FILES ${filesPYTests} DESTINATION ${SUDODEM_PY_PATH}/sudodem/tests)
#FILE(GLOB filesPYPerf "${CMAKE_CURRENT_SOURCE_DIR}/../examples/test/performance/*")
#INSTALL(FILES ${filesPYPerf} DESTINATION ${SUDODEM_PY_PATH}/sudodem/tests/checks/performance)

ADD_LIBRARY(WeightedAverage2d SHARED "${CMAKE_CURRENT_SOURCE_DIR}/WeightedAverage2d.cpp")
SET_TARGET_PROPERTIES(WeightedAverage2d PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../")
TARGET_LINK_LIBRARIES(WeightedAverage2d sudodem)
INSTALL(TARGETS WeightedAverage2d DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

#FILE(GLOB filesFEMxDEM "${CMAKE_CURRENT_SOURCE_DIR}/FEMxDEM/*.py")
#INSTALL(FILES ${filesFEMxDEM} DESTINATION ${SUDODEM_PY_PATH}/sudodem/FEMxDEM)

ADD_LIBRARY(_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_utils.cpp")
SET_TARGET_PROPERTIES(_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../")
TARGET_LINK_LIBRARIES(_utils sudodem)
INSTALL(TARGETS _utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")


ADD_LIBRARY(_polyhedra_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_polyhedra_utils.cpp")
SET_TARGET_PROPERTIES(_polyhedra_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../")
TARGET_LINK_LIBRARIES(_polyhedra_utils sudodem)
INSTALL(TARGETS _polyhedra_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")


#ADD_LIBRARY(_packPredicates SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packPredicates.cpp")
#SET_TARGET_PROPERTIES(_packPredicates PROPERTIES PREFIX "")
#TARGET_LINK_LIBRARIES(_packPredicates sudodem)
#IF(ENABLE_GTS AND NOT(PY_gts))
#  TARGET_LINK_LIBRARIES(_packPredicates _gts)
#ENDIF(ENABLE_GTS AND NOT(PY_gts))


#INSTALL(TARGETS _packPredicates DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

#ADD_LIBRARY(_packDisks SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packDisks.cpp")
#SET_TARGET_PROPERTIES(_packDisks PROPERTIES PREFIX "")
#TARGET_LINK_LIBRARIES(_packDisks sudodem)
#INSTALL(TARGETS _packDisks DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

#ADD_LIBRARY(_packObb SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packObb.cpp")
#SET_TARGET_PROPERTIES(_packObb PROPERTIES PREFIX "")
#INSTALL(TARGETS _packObb DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(wrapper SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/sudodemWrapper.cpp")
SET_TARGET_PROPERTIES(wrapper PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(wrapper sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS wrapper DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_customConverters SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/customConverters.cpp")
SET_TARGET_PROPERTIES(_customConverters PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_customConverters sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _customConverters DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_superellipse_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_superellipse_utils.cpp;${CMAKE_CURRENT_SOURCE_DIR}/svgobjects.cpp")
SET_TARGET_PROPERTIES(_superellipse_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_superellipse_utils sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _superellipse_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

#ADD_LIBRARY(_fem_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_fem_utils.cpp")
#SET_TARGET_PROPERTIES(_fem_utils PROPERTIES PREFIX "")
#TARGET_LINK_LIBRARIES(_fem_utils)
#INSTALL(TARGETS _fem_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

#ADD_LIBRARY(_gjkparticle_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_gjkparticle_utils.cpp")
#SET_TARGET_PROPERTIES(_gjkparticle_utils PROPERTIES PREFIX "")
#TARGET_LINK_LIBRARIES(_gjkparticle_utils voro++)
#INSTALL(TARGETS _gjkparticle_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")
