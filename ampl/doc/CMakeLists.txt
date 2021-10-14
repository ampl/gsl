# Check requirements.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/support) # for FindSphinx
find_package(Python)
find_package(Sphinx)
if(NOT SPHINX_EXECUTABLE)
  message(STATUS "Documentation generation disabled (requires sphinx)")
endif()

set(BUILD_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/support/extract-docs.py)

# Get .rst file names from amplgsl.cc.
set(docsStaticDir ${CMAKE_CURRENT_SOURCE_DIR}/amplgsl)
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

# Add static documents to amplgsl_docs to then add them as dependecies to the target
# and ensure they get built
set(staticDocs accuracy examples fdl  freedoc front-matter  gpl  history 
      index  intro no-warranty  rng)
foreach(f ${staticDocs})
    set(amplgsl_docs ${amplgsl_docs} ${docsStaticDir}/${f}.rst)
endforeach()

set(DOCSBASEDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_SOURCE  ${CMAKE_CURRENT_SOURCE_DIR}/amplgsl)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/output)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/xml/index.xml)

file(GLOB_RECURSE DOC_SOURCES ${SPHINX_SOURCE}/*.rst)

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
                  ${DOCSBASEDIR} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS ${amplgsl_docs}
                  MAIN_DEPENDENCY ${DOCSBASEDIR}/conf.py
                  COMMENT "Generating documentation with Sphinx")

# Nice named target so we can run the job easily
add_custom_target(sphinx ALL DEPENDS ${SPHINX_INDEX_FILE}
SOURCES ${DOC_SOURCES} )
add_to_folder(docs sphinx doxygen)




#add_custom_target(amplgsl-doc
  #COMMAND ${Python_EXECUTABLE} ${BUILD_DOCS}
  #DEPENDS conf.py ${amplgsl_docs})

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ampl.github.io/
        DESTINATION amplgsl OPTIONAL
        PATTERN .* EXCLUDE)
