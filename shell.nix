with (import <nixpkgs> {});



let
  mip = ps: ps.callPackage ./mip.nix {};

  pythonEnv = python38.withPackages (
      ps: [ ps.pymatgen
            ps.monty
            ps.openbabel-bindings
            ps.pygraphviz
            ps.mpi4py
            (mip ps)
          ]);
in mkShell rec {
  buildInputs = [ pythonEnv
                  sqlitebrowser
                  texlive.combined.scheme-small
                  mpi
                  cbc
                  (import ./RNMC.nix)
                ];
}
