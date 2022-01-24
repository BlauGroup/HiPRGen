{
  description = "Multi node reaction network generator";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-21.05;
    RNMC.url = github:BlauGroup/RNMC;
    flake-compat = {
      url = github:edolstra/flake-compat;
      flake = false;
    };
  };

  outputs = { self, nixpkgs, RNMC, flake-compat }:

    let
      HiPRGen = systemString:
        with import nixpkgs { system = systemString; };
        with python38Packages;
        buildPythonPackage {
          pname = "HiPRGen";
          version = "0.2";
          src = ./.;
          checkInputs = [
            pymatgen
            monty
            openbabel-bindings
            pygraphviz
            mpi4py
            pycairo
            mpi
            (builtins.getAttr systemString RNMC.defaultPackage)
            sqlite
            openssh # needed for correct MPI functioning
          ];

          checkPhase = "python test.py 2";
        };


      genericDevShell = systemString: installHiPRGen:
        with import nixpkgs { system = systemString; };
        mkShell {
          buildInputs = [
            (python38.withPackages (
              ps: [ ps.pymatgen
                    ps.monty
                    ps.openbabel-bindings
                    ps.pygraphviz
                    ps.mpi4py
                    ps.pycairo
                    (if installHiPRGen then (HiPRGen systemString) else null)
                    ps.jupyterlab
                  ]))

            texlive.combined.scheme-small
            mpi
            (sqlite.override { interactive = true; })
            (builtins.getAttr systemString RNMC.defaultPackage)
          ];
        };

    in {
      devShell = {
        x86_64-linux = genericDevShell "x86_64-linux" false;
        x86_64-darwin = genericDevShell "x86_64-darwin" false;
        aarch64-darwin = genericDevShell "aarch64-darwin" false;
      };

      defaultPackage = {
        x86_64-linux = HiPRGen "x86_64-linux";
        x86_64-darwin = HiPRGen "x86_64-darwin";
        aarch64-darwin = HiPRGen "aarch64-darwin";
      };

      checks = {
        x86_64-linux.tests = HiPRGen "x86_64-linux";
        x86_64-darwin.tests = HiPRGen "x86_64-darwin";
        aarch64-darwin = HiPRGen "aarch64-darwin";
      };
    };

}
