function(generate_clangd_file COMPILER_FLAGS)

  set(CLANGD_COMPILER_STANDARD "-std=c++17")

  string(CONCAT CLANGD_FILE_CONTENT 
    "CompileFlags:\n"
    "  CompilationDatabase: ${CMAKE_BINARY_DIR}\n"
    "  Add:\n"
    "    - "
    ${CLANGD_COMPILER_STANDARD}
    "\n")
  
    foreach(FLAG ${COMPILER_FLAGS})
      string(CONCAT CLANGD_FILE_CONTENT 
      "${CLANGD_FILE_CONTENT}"
      "    - "
      ${FLAG}
      "\n")
    endforeach()
    

  file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/.clangd "${CLANGD_FILE_CONTENT}")
endfunction()
