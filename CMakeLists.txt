project(bcalm)

cmake_minimum_required(VERSION 2.6)

################################################################################
# Shortcuts
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/gatb-core)

################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${GATB_CORE_HOME}/gatb-core/cmake")
    IF (EXISTS "${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${path}")
    ENDIF()
ENDFOREACH(path)

#############################
#getting git version
#from http://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake
exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "rev-parse --short HEAD"
    OUTPUT_VARIABLE VERSION_SHA1
    RETURN_VALUE ERROR_GIT)

if(NOT ${ERROR_GIT})
    add_definitions( -DGIT_SHA1="${VERSION_SHA1}" )
else()
    message("Warning: cannot retrieve git version. Bcalm won't display its version. Error value: ${ERROR_GIT})")
endif(NOT ${ERROR_GIT})

################################
#add version nifo

file (STRINGS "VERSION" VERSION)
add_definitions( -DVERSION="${VERSION}" )

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1) # no need to compile dbgh5, etc..
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)


# GATB CORE
include (GatbCore)

################################################################################
#  TOOLS
################################################################################

MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET (dirlist "")
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
          LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    SET(${result} ${dirlist})
ENDMACRO()

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

set (program "bcalm")
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
add_executable(${program} ${ProjectFiles})
target_link_libraries(${program} ${gatb-core-libraries})

################################################################################
# DELIVERY
################################################################################

# If your current login name is different from your GForge login name, you have
# to overwrite the CPACK_USER_NAME to be the same as your GForge login
#SET (CPACK_USER_NAME                    "your_gforge_login")

# We set the version number
SET (CPACK_PACKAGE_VERSION              "")  

# We have to tell what is the server name
SET (CPACK_GFORGE_PROJECT_NAME          "gatb-tools")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwanted files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES
    "^${PROJECT_SOURCE_DIR}/.git/"     ;
    "^${PROJECT_SOURCE_DIR}/.gitmodules" ;
    "^${PROJECT_SOURCE_DIR}/.gitignore" ;
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/.cproject" ;
    "^${GATB_CORE_HOME}/.git/"     ;
    "^${GATB_CORE_HOME}/.project"  ;
    "^${GATB_CORE_HOME}/.gitignore"
)

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
INSTALL (DIRECTORY "${PROJECT_SOURCE_DIR}/example/" DESTINATION example)
INSTALL (FILES VERSION LICENSE README.md DESTINATION bin/..)

# cmake_system_name for mac is Darwin, i want to change that
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(PRETTY_SYSTEM_NAME "Linux")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(PRETTY_SYSTEM_NAME "Mac")
endif()

set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-binaries-${VERSION}-${PRETTY_SYSTEM_NAME})

include (CPack)


