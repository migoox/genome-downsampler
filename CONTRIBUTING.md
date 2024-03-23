# Contributing

## File structure

```
Foo
 ├── include
 │   └── Foo
 │       ├── Foo.h
 │       └── ...
 ├── src
 │   ├── Foo.cpp
 │   └── ...
 ├── test
 │   └── ...
 └── libs
     ├── A
     │   ├── include
     │   │   └── A
     │   │       ├── A.h
     │   │       └── ...
     │   ├── src
     │   │   ├── A.cpp
     │   │   └── ...
     │   └── test
     │       └── ...
     └── B
         ├── include
         │   └── B
         │       ├── B.h
         │       └── ...
         ├── src
         │   ├── B.cpp
         │   └── ...
         └── test
             └── ...
```

https://github.com/Jamagas/CMake

## Setup

Required programs:

- clang/clang++ (at least 14.0) – LLVM’s frontend for g++ and msvc (it overlays both), which helps us with cross-platform apps,
- clangd (at least 14.0) – LLVM’s c/cpp language server, allows auto-completion and clang-tidy features,
- clang-format (at least 14.0) – LLVM’s c/cpp code formatter,
- stdlibc++ (for clang 14.0 we need at least 12.0) – GNU implementation of C++ standard library, there is also LLVM’s libc++, but install stdlibc++,
- cmake (at least 3.21),
- cmake-format (at least 0.6) – allows cmake files formatting,

Required vscode extensions:

- clang-format (author: Xaver Hellauer) – provides clang formatting integration,
- clangd (author: LLVM) – clangd is a language server, and provides C++ IDE features to editors,
- cmake (author: twxs) – provides autocompletion, and syntax highliting for cmake,
- cmake tools (author: Microsoft) – provides cmake commands for build, run, test etc.
- cmake-format (author: cheshirekow) – provides cmake formatting integration,

VSCode configuration:

```json
// Autoformatting
  "editor.formatOnSave": true,

  "[cpp]": {
// Use clang formatter
    "editor.defaultFormatter": "xaver.clang-format"
  },

  "cmake.options.statusBarVisibility": "visible",

// Set clangd as default provider
  "C_Cpp.default.configurationProvider": "clangd",
// Disable Microsoft's intelli sense
  "C_Cpp.intelliSenseEngine": "disabled",

// Find your executables for these programs (which command)
  "clang-format.executable": "/usr/bin/clang-format",
  "clangd.path": "/usr/bin/clangd",
  "cmakeFormat.exePath": "/usr/bin/cmake-format",

  "[cmake]": {
// Use cmake formatter
    "editor.defaultFormatter": "cheshirekow.cmake-format"
  },
```

VSCode how to use:
When the project is opened for the first time vscode will ask you for choosing cmake kit – check clang kit.

Don’t use regular run, since it will use g++, when you want to build or run your app use `Ctrl+Shfig+P` and choose CMake Build/Run/Debug.