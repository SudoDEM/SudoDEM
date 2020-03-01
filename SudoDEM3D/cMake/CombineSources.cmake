# include several source files (in SRCS) in one or more
# combined files; each combined file holds maximally 
# MAXNUM files.
# BASE gives basename for created files;
# files corresponding to the pattern are deleted
# first, so that there are no relicts from previous runs
# with possibly different MAXNUM
MACRO(COMBINE_SOURCES BASE SRCS MAXNUM)
	LIST(LENGTH SRCS SRCS_LENGTH)
	SET(COMB_COUNTER 0)
	FILE(GLOB EXISTING "${BASE}.*.cpp")
	IF("$EXISTING")
		FILE(REMOVE ${EXISTING})
	ENDIF()
	SET(OUT "${BASE}.${COMB_COUNTER}.cpp")
	FILE(WRITE ${OUT})
	SET(COUNTER 0)
	FOREACH(SRC ${SRCS})
		if(${SRC} MATCHES "^/.*$") # absolute filename
			FILE(APPEND ${OUT} "#include<${SRC}>\n")
		else()
			FILE(APPEND ${OUT} "#include<${CMAKE_SOURCE_DIR}/${SRC}>\n")
		endif()
		MATH(EXPR COUNTER "${COUNTER}+1")
		IF(${COUNTER} EQUAL ${MAXNUM})
			SET(COUNTER 0)
			MATH(EXPR COMB_COUNTER ${COMB_COUNTER}+1)
			SET(OUT "${BASE}.${COMB_COUNTER}.cpp")
			FILE(WRITE ${OUT})
		ENDIF()
	ENDFOREACH()
ENDMACRO()
