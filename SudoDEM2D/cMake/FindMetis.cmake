# - Find Metis library
# 
# This module defines
#  METIS_INCLUDE_DIR, where to find loki/Typelist.h, etc.
#  METIS_LIBRARY, libraries to link against to use GL2PS.
#  METIS_FOUND, If false, do not try to use GL2PS.

FIND_PATH(METIS_INCLUDE_DIR metis.h PATHS /usr/include/metis)
FIND_LIBRARY(METIS_LIBRARY NAMES metis parmetis PATHS /usr/lib)

# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Metis  DEFAULT_MSG  METIS_INCLUDE_DIR METIS_LIBRARY)

MARK_AS_ADVANCED(METIS_INCLUDE_DIR METIS_LIBRARY)
