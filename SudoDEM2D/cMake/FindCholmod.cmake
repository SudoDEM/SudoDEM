# - Try to find CHOLMOD
# This will define
#
#  CHOLMOD_FOUND          - system has CHOLMOD
#  CHOLMOD_LIBRARIES 	    - library to link against to use Cholmod 
#  CHOLMOD_INCLUDE_DIR    - where to find cholmod.h, etc.
#  AMD_LIBRARY	 	        - needed by CHOLMOD
#  COLAMD_LIBRARY 	      - needed by CHOLMOD
#  CCOLAMD_LIBRARY 	      - needed by CHOLMOD
#  CAMD_LIBRARY 	        - needed by CHOLMOD

FIND_LIBRARY(CHOLMOD_LIBRARIES NAMES cholmod libcholmod
        PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib/CGAL
        /usr/lib64
        /usr/local/lib64
        /usr/lib64/CGAL
    )

FIND_LIBRARY(AMD_LIBRARY NAMES amd PATHS /usr/lib /usr/local/lib /usr/lib/CGAL /usr/lib64 /usr/local/lib64 /usr/lib64/CGAL)
FIND_LIBRARY(CAMD_LIBRARY NAMES camd PATHS /usr/lib /usr/local/lib /usr/lib/CGAL /usr/lib64 /usr/local/lib64 /usr/lib64/CGAL)
FIND_LIBRARY(COLAMD_LIBRARY NAMES colamd PATHS /usr/lib /usr/local/lib /usr/lib/CGAL /usr/lib64 /usr/local/lib64 /usr/lib64/CGAL)
FIND_LIBRARY(CCOLAMD_LIBRARY NAMES ccolamd PATHS /usr/lib /usr/local/lib /usr/lib/CGAL /usr/lib64 /usr/local/lib64 /usr/lib64/CGAL)

FIND_PATH(CHOLMOD_INCLUDE_DIR cholmod.h PATH /usr/include /usr/include/suitesparse)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cholmod DEFAULT_MSG CHOLMOD_LIBRARIES CHOLMOD_INCLUDE_DIR AMD_LIBRARY CAMD_LIBRARY COLAMD_LIBRARY CCOLAMD_LIBRARY)
MARK_AS_ADVANCED(CHOLMOD_LIBRARIES CHOLMOD_INCLUDE_DIR AMD_LIBRARY CAMD_LIBRARY COLAMD_LIBRARY CCOLAMD_LIBRARY)
