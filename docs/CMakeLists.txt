find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# Find all the public headers
get_target_property(QSF_PUBLIC_HEADER_DIR QSF INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE QSF_PUBLIC_HEADERS ${QSF_PUBLIC_HEADER_DIR}/*.h)

set(DOXYGEN_INDEX_FILE xml/index.xml)

# file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${QSF_PUBLIC_HEADERS}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   COMMAND ${DOXYGEN_EXECUTABLE} 
                   COMMENT "Generating documentation with Doxygen")

add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/rst)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)
# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   DEPENDS ${SPHINX_SOURCE}/index.rst ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx")
# Nice named target so we can run the job easily
add_custom_target(sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})
