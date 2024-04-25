{
  description = "Nix project environment with clang, cmake and tooling";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
  };

  outputs = { self, nixpkgs }:
    let
      supportedSystems = [ "x86_64-linux" ];
      forAllSystems = f: nixpkgs.lib.genAttrs supportedSystems (system: f system);
    in
    {
      devShells = forAllSystems
        (system:
          let
            pkgs = import nixpkgs {
              inherit system;
              config.allowUnfree = true;
            };
          in
          {
            default = with pkgs; mkShell.override { stdenv = gcc12Stdenv; }
              {
                shellHook = ''
                  export CUDA_PATH=${cudatoolkit}
                  export LD_LIBRARY_PATH=${linuxPackages.nvidia_x11}/lib:${ncurses5}/lib
                  export EXTRA_LDFLAGS="-L/lib -L${linuxPackages.nvidia_x11}/lib"
                  export EXTRA_CCFLAGS="-I/usr/include"
                  export NVCC_APPEND_FLAGS="-L${gcc13.cc.lib}/lib"
                  export HTLSIB_ROOT=${htslib}

                  zellij
                '';
                packages = [
                  clang-tools
                  bear
                  cmake
                  cmake-format
                  zellij
                  pkg-config
                ];
                buildInputs = [
                  cudatoolkit linuxPackages.nvidia_x11
                  ncurses5 binutils

                  boost htslib

                  # or-tools
                  or-tools
                  bzip2 cbc eigen glpk python3.pkgs.absl-py
                  python3.pkgs.pybind11 python3.pkgs.setuptools
                  python3.pkgs.wheel re2 zlib
                ];
              };
          }
        );
    };
}
