with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "d4eb6aeca344ec4fb775a68e8dd6d7bb03a43951";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "1pzldq21pqs545ykn4jraaqrs2lv9vr1706vym0g9fa7bj99da09";
  };
  buildPhase = "CC=gcc ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}
