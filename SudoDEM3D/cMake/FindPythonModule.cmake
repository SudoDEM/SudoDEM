# http://www.cmake.org/pipermail/cmake/2011-January/041666.html
#
# - Find Python Module

FUNCTION(find_python_module module)
  STRING(TOUPPER ${module} module_upper)
  
  IF(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
    SET(${module}_FIND_REQUIRED TRUE)
  ENDIF(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
  
  EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
    "import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
    RESULT_VARIABLE _${module}_status 
    OUTPUT_VARIABLE _${module}_location
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    
    IF(_${module}_status MATCHES 0)
      SET(PY_${module} ${_${module}_location} CACHE STRING "Location of Python module ${module}")
    ENDIF(_${module}_status MATCHES 0)
    
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(${module} DEFAULT_MSG PY_${module})
ENDFUNCTION(find_python_module)
