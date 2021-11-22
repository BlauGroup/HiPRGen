{
  description = "Multi node reaction network generator";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-21.05;
    RNMC.url = github:BlauGroup/RNMC;
  };

  outputs = { self, nixpkgs }:
    let genericDevShell = systemString:
          with import nixpkgs { system = systemString; };
          let
            HiPRGen = python38Packages.buildPythonPackage {
              pname = "HiPRGen";
              version = "0.1";
              src = ./.;
              doCheck = false;
            };

            pythonEnv = python38.withPackages (
              ps: [ ps.pymatgen
                    ps.monty
                    ps.openbabel-bindings
                    ps.pygraphviz
                    ps.mpi4py
                    ps.pycairo
                    # HiPRGen
                  ]);

          in mkShell {
            buildInputs = [ pythonEnv
                            texlive.combined.scheme-small
                            mpi
                            sqlite
                            RNMC
                          ];
          };
    in
      devShell = {
        x86_64-linux = genericDevShell "x86_64-linux";
      };

}
