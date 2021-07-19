with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "1f948bc009696f9b06d27fb60b994b1e94f1edc0";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "15d67s1g2aj5ln4gcbd45lv5s88djgi5rrkm0bvslv9dsdjfbxjb";
  };
  buildPhase = "CC=gcc ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}


