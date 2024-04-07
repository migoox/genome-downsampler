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
                  export LD_LIBRARY_PATH=${linuxPackages.nvidia_x11}/lib:${pkgs.ncurses5}/lib
                  export EXTRA_LDFLAGS="-L/lib -L${linuxPackages.nvidia_x11}/lib"
                  export EXTRA_CCFLAGS="-I/usr/include"
                  export HTLSIB_ROOT=${htslib}

                  zellij
                '';
                packages = [
                  clang-tools
                  zellij
                  bear
                  cmake
                  cmake-format
                  zellij
                ];
                buildInputs = [
                  llvmPackages.libcxxClang

                  cudatoolkit linuxPackages.nvidia_x11
                  ncurses5 binutils

                  boost htslib
                ];
              };
          }
        );
    };
}
