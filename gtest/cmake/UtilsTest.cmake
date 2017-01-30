##---------------------------------------------------------------------------##
## src/cmake/UtilsTest.cmake
## Thomas M. Evans
##---------------------------------------------------------------------------##
## Setup Utils test expressions and services
##---------------------------------------------------------------------------##
## Setup PASS/FAIL expressions for Utils test harness

# Utils-test-harness pass/fail
SET(UtilsPass "Test: PASSED")
SET(UtilsFail "Test: FAILED")

# Utils google test harness pass/fail
SET(UtilsGtestPass "overall test result: PASSED")
SET(UtilsGtestFail "overall test result: FAILED")

# Turn off variadic-macros warnings when using gtest
IF (CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
ENDIF()

# Set number of parallel tests to run
IF (TPL_ENABLE_MPI)
  SET(UtilsNP 1 2 4)
ELSE()
  SET(UtilsNP "1")
ENDIF()

INCLUDE(CMakeParseArguments)

##---------------------------------------------------------------------------##
## ADDING UNIT TEST
##---------------------------------------------------------------------------##
# ADD_UTILS_TEST(
#   SOURCE_FILE
#   [NP 1 [2 [...]]]
#   [DEPLIBS lib1 [lib2 ...]]
#   [ENVIRONMENT VAR=value [VAR2=value2 ...]]
#   [ISOLATE]
#   [LEGACY]
#   [DISABLE]
#   )
#
# Create and add a unit test from the source file SOURCE_FILE.
#
# NP specifies the number of processors to use for this unit test. The default
# is to use UtilsNP (1, 2, and 4) for MPI builds and 1 for serial builds.
#
# DEPLIBS specifies extra libraries to link to. By default, unit tests will link
# against the package's current library; Gtest unit tests will also link against
# Utils_gtest.
#
# ENVRIONMENT sets the given environmental variables when the test is run.
#
# If ISOLATE is specified, the test will be run in its own directory.
#
# If LEGACY is specified, we use the old Utils harness. Otherwise, use the
# Utils google test harness.
#
# If DISABLE is specified, we will build the test executable but omit it from
# the list of tests to run through CTest.
#
FUNCTION(ADD_UTILS_TEST SOURCE_FILE)
  cmake_parse_arguments(PARSE
    "LEGACY;ISOLATE;DISABLE"
    ""
    "DEPLIBS;NP;ENVIRONMENT" ${ARGN})

  # Set googletest/harness options
  IF (PARSE_LEGACY)
    MESSAGE(ERROR "Legacy testing unavailable in ExaMPM.")
  ELSE()
    SET(PASS_RE ${UtilsGtestPass})
    SET(FAIL_RE ${UtilsGtestFail})
  ENDIF()

  # Add additional library dependencies if needed
  IF(PARSE_DEPLIBS)
    SET(DEPLIBS_CMD TESTONLYLIBS ${PARSE_DEPLIBS})
  ENDIF()

  # Set number of processors, defaulting to UtilsNP
  SET(NUM_PROCS ${PARSE_NP})
  IF (NOT NUM_PROCS)
    SET(NUM_PROCS ${UtilsNP})
  ENDIF()

  # Check to see if MPI-only unit test
  LIST(FIND NUM_PROCS 1 HAS_SERIAL)
  IF (HAS_SERIAL EQUAL -1)
    SET(COMM mpi)
    IF(NOT TPL_ENABLE_MPI)
      # return early to avoid potential set_property on nonexistent test
      RETURN()
    ENDIF()
  ELSE()
    SET(COMM serial mpi)
  ENDIF()

  # add the test executable
  GET_FILENAME_COMPONENT(EXE_NAME ${SOURCE_FILE} NAME_WE)
  TRIBITS_ADD_EXECUTABLE(
    ${EXE_NAME}
    SOURCES ${SOURCE_FILE}
    ${DEPLIBS_CMD}
    COMM ${COMM}
    )

  # If the test is disabled, print a small message and omit it from the CTest
  IF (PARSE_DISABLE)
    MESSAGE("Disabling testing for ${SOURCE_FILE} in ${SUBPACKAGE_FULLNAME}")
    RETURN()
  ENDIF()

  # Loop over processors for parallel tests
  FOREACH(np ${NUM_PROCS})
    IF (PARSE_ISOLATE)
      # Add an "advanced" test
      TRIBITS_ADD_ADVANCED_TEST(
        ${EXE_NAME}_MPI_${np}
        TEST_0
          EXEC ${EXE_NAME}
          NUM_MPI_PROCS ${np}
          PASS_REGULAR_EXPRESSION "${PASS_RE}"
          FAIL_REGULAR_EXPRESSION "${FAIL_RE}"
        OVERALL_WORKING_DIRECTORY TEST_NAME
        )
    ELSE()
      # Add a normal test
      TRIBITS_ADD_TEST(
        ${EXE_NAME}
        NUM_MPI_PROCS ${np}
        PASS_REGULAR_EXPRESSION "${PASS_RE}"
        FAIL_REGULAR_EXPRESSION "${FAIL_RE}"
        )
    ENDIF()
  ENDFOREACH()

  # set environmental variables if necessary
  IF(PARSE_ENVIRONMENT)
    FOREACH(np ${NUM_PROCS})
      IF (TPL_ENABLE_MPI)
        SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${EXE_NAME}_MPI_${np}")
      ELSE()
        SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${EXE_NAME}")
      ENDIF()

      # Modify environment
      SET_PROPERTY(TEST "${TEST_NAME}"
        PROPERTY ENVIRONMENT
        ${PARSE_ENVIRONMENT}
        )
    ENDFOREACH()
  ENDIF()

ENDFUNCTION()

##---------------------------------------------------------------------------##
##                   end of UtilsTest.cmake
##---------------------------------------------------------------------------##
