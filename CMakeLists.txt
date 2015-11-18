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
# DELIVERY
################################################################################

# If your current login name is different from your GForge login name, you have
# to overwrite the CPACK_USER_NAME to be the same as your GForge login
#SET (CPACK_USER_NAME                    "your_gforge_login")

# We set the version number
SET (CPACK_PACKAGE_VERSION              "1.0.0")  

# We have to tell what is the server name
SET (CPACK_GFORGE_PROJECT_NAME          "gatb-tools")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwated files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
    "^${GATB_CORE_HOME}/.project"
    "^${GATB_CORE_HOME}/.gitignore"
    "^${GATB_CORE_HOME}/doc/design"
    "^${GATB_CORE_HOME}/DELIVERY.md"
)

include (Delivery)


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

set (program "bglue")
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/bglue)
include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles    FOLLOW_SYMLINKS ${PROGRAM_SOURCE_DIR}/*.cpp )
add_executable(${program} ${ProjectFiles})
target_link_libraries(${program} ${gatb-core-libraries})


