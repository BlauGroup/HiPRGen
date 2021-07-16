with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "b79436dc27feff194e485de649131784d7a6931c";
  buildInputs = [gsl clang sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "17md2pjsxlwg3cd5blr2fzjr0x369z9xsza4imyn7g6l67x8ym9g";
  };
  buildPhase = "CC=clang ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}

