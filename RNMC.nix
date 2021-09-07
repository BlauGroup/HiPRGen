with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "eac44cc9d8edd79e3cda9b29290341d6df447f2c";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "0jr60yzcz3ccmgkc3cl6bckpskf271wjn310nz4rydpxw4941yw3";
  };
  buildPhase = "CC=gcc ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}
