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
            };
          in
          {
            default = with pkgs; mkShell.override { stdenv = clangStdenv; }
              {
                shellHook = ''
                  zellij
                '';
                packages = [
                  clang-tools
                  zellij
                  bear
                  cmake
                  cmake-format
                ];
                buildInputs = [
                  llvmPackages.libcxxClang
                ];
              };
          }
        );
    };
}
