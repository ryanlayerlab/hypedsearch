set(CRUX_VERSION_FULL "${CRUX_VERSION}")

if (EXISTS "${SOURCE_DIR}/.git/HEAD")
  # Add the git hash
  execute_process(
    COMMAND "git" "rev-parse" "--short" "HEAD"
    WORKING_DIRECTORY ${SOURCE_DIR}
    RESULT_VARIABLE GIT_RESULT
    OUTPUT_VARIABLE PROJECT_SOURCE_VERSION
  )
  if (GIT_RESULT EQUAL 0)
    string(
      REGEX REPLACE "\n$" ""
      PROJECT_SOURCE_VERSION
      ${PROJECT_SOURCE_VERSION}
    )
    set(CRUX_VERSION_FULL "${CRUX_VERSION_FULL}-${PROJECT_SOURCE_VERSION}")
  else (GIT_RESULT EQUAL 0)
      MESSAGE("Git command git show-ref failed: ${PROJECT_SOURCE_VERSION}")
  endif (GIT_RESULT EQUAL 0)

  # Add the last commit date
  execute_process(
    COMMAND "git" "log" "-1" "--format=%cd" "--date=short"
    WORKING_DIRECTORY ${SOURCE_DIR}
    RESULT_VARIABLE GIT_DATE_RESULT
    OUTPUT_VARIABLE PROJECT_COMMIT_DATE
  )
  if (GIT_DATE_RESULT EQUAL 0)
    string(
      REGEX REPLACE "\n$" ""
      PROJECT_COMMIT_DATE
      ${PROJECT_COMMIT_DATE}
    )
    set(CRUX_VERSION_FULL "${CRUX_VERSION_FULL}-${PROJECT_COMMIT_DATE}")
  else (GIT_DATE_RESULT EQUAL 0)
    MESSAGE("Git command git show-ref failed: ${PROJECT_COMMIT_DATE}")
  endif (GIT_DATE_RESULT EQUAL 0)

endif (EXISTS "${SOURCE_DIR}/.git/HEAD")


file(WRITE crux_version.h.txt "#define CRUX_VERSION \"${CRUX_VERSION_FULL}\"\n")

# Copy the file to the final header if the version changed
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different crux_version.h.txt crux_version.h)

