# - Find OpenBlas library
# 
# This module defines
#  OPENBLAS_LIBRARY, libraries to link against to use Openblas.
#  OPENBLAS_FOUND, If false, do not try to use Openblas.

FIND_LIBRARY(OPENBLAS_LIBRARY NAMES openblas blas PATHS /usr/lib/openblas-base )

# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenBlas  DEFAULT_MSG  OPENBLAS_LIBRARY)

MARK_AS_ADVANCED(OPENBLAS_LIBRARY)
