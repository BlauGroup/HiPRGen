with (import <nixpkgs> {});


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
                  (import ./RNMC.nix)
                ];
}
