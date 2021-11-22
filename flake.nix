{
  description = "Multi node reaction network generator";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-21.05;
    RNMC.url = github:BlauGroup/RNMC;
  };

  outputs = { self, nixpkgs, RNMC }:


    let genericDevShell = systemString: includeSelf:
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
                    (if includeSelf then HiPRGen else null)
                  ]);

          in mkShell {
            buildInputs = [ pythonEnv
                            texlive.combined.scheme-small
                            mpi
                            sqlite
                            (builtins.getAttr systemString RNMC.defaultPackage)
                          ];
          };
    in {
      devShell = {
        x86_64-linux = genericDevShell "x86_64-linux" false;
      };
    };

}
