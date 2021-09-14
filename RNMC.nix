with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "869021cc4ca19757ac3b6617b95c7ee6c81b6cf3";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "0n7mw7w95ypvi7gxf65kv2r8r2z4c2w2mahfw2r5sj0f890b8dhi";
  };
  buildPhase = "CC=gcc ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}
