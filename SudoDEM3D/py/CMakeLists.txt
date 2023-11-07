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

#==================pyModules=========================================

FILE(GLOB filesPY "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
INSTALL(FILES ${filesPY} DESTINATION ${SUDODEM_PY_PATH}/sudodem)
INSTALL(FILES 3rd-party/mtTkinter-0.4/mtTkinter.py DESTINATION ${SUDODEM_LIB_PATH}/py)

ADD_LIBRARY(_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_utils.cpp")
SET_TARGET_PROPERTIES(_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/qt;$ORIGIN/../../")
TARGET_LINK_LIBRARIES(_utils sudodem _GLViewer)
INSTALL(TARGETS _utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(wrapper SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/sudodemWrapper.cpp")
SET_TARGET_PROPERTIES(wrapper PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(wrapper sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS wrapper DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_customConverters SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/customConverters.cpp")
SET_TARGET_PROPERTIES(_customConverters PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_customConverters sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _customConverters DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_superquadrics_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_superquadrics_utils.cpp")
SET_TARGET_PROPERTIES(_superquadrics_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_superquadrics_utils sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _superquadrics_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_fem_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_fem_utils.cpp")
SET_TARGET_PROPERTIES(_fem_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_fem_utils sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _fem_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")

ADD_LIBRARY(_gjkparticle_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_gjkparticle_utils.cpp")
SET_TARGET_PROPERTIES(_gjkparticle_utils PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;$ORIGIN/../../;$ORIGIN/../../../3rdlibs")
TARGET_LINK_LIBRARIES(_gjkparticle_utils voro++ sudodem ${Boost_LIBRARIES})
INSTALL(TARGETS _gjkparticle_utils DESTINATION "${SUDODEM_PY_PATH}/sudodem/")
