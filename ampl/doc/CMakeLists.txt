# Check requirements.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/source/support) # for FindSphinx
find_package(Python)
find_package(Sphinx)
if(NOT SPHINX_EXECUTABLE)
  message(STATUS "Documentation generation disabled (requires sphinx)")
  return()
endif()


set(BUILD_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/source/support/extractdocs.py)

# Get .rst file names from amplgsl.cc.
set(docsStaticDir ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(docsDestinationDir ${docsStaticDir}/ref)

set(amplgsl_docs )
set(amplgsl_src ${CMAKE_CURRENT_SOURCE_DIR}/../src/amplgsl.cc)
file(READ ${amplgsl_src} content)
string(REGEX MATCHALL "@file[^\n]*" files ${content})
foreach (f ${files})
   if (f MATCHES "@file +(.*)")
     set(amplgsl_docs ${amplgsl_docs}
         ${docsDestinationDir}/${CMAKE_MATCH_1}.rst)
   endif ()
endforeach ()
# Add a command to extract the amplgsl documentation from amplgsl.cc.
add_custom_command(OUTPUT ${amplgsl_docs}
  COMMAND ${Python_EXECUTABLE} ${BUILD_DOCS} ${amplgsl_src} -o ${docsDestinationDir}
  DEPENDS ${BUILD_DOCS} ${amplgsl_src})

set(DOCSBASEDIR ${CMAKE_CURRENT_SOURCE_DIR}/source/)
set(SPHINX_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/output)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/html/index.html)

file(GLOB DOC_SOURCES ${DOCSBASEDIR}/*.rst)

add_to_source_group(REF_SOURCES ref "" ${amplgsl_docs})
add_to_source_group(DOC_SOURCES docs "" ${DOC_SOURCES})
set(DOC_INDEX ${DOCSBASEDIR}/index.rst)

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
                  ${DOCSBASEDIR} ${SPHINX_OUTPUT_DIR}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS ${amplgsl_docs} ${DOC_SOURCES} ${DOC_INDEX}
                  MAIN_DEPENDENCY ${DOCSBASEDIR}/conf.py
                  COMMENT "Generating documentation with Sphinx")

# Nice named target so we can run the job easily
add_custom_target(doc ALL DEPENDS ${SPHINX_INDEX_FILE}
SOURCES ${REF_SOURCES} ${DOC_SOURCES}  ${DOC_INDEX})

add_to_folder(docs doc)


install(DIRECTORY ${SPHINX_OUTPUT_DIR}/
        DESTINATION amplgsl/doc OPTIONAL
        PATTERN .* EXCLUDE)
