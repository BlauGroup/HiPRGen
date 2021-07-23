with (import <nixpkgs> {});
let
  python = let
    packageOverrides = self: super: {
      pyarrow = super.pyarrow.overridePythonAttrs (
        old: { PYARROW_WITH_PLASMA = true;}
      );


    }; in python38.override {inherit packageOverrides;};

  pythonEnv = python.withPackages (
      ps: [ ps.pymatgen
            ps.monty
            ps.openbabel-bindings
            ps.pygraphviz
            ps.mpi4py
            ps.pyarrow
          ]);
in mkShell rec {
  buildInputs = [ pythonEnv
                  texlive.combined.scheme-small
                  mpi
                  (import ./RNMC.nix)
                ];
}
