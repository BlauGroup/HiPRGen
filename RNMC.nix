with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "c770edb9606f07c26438e61936a375a686c2994c";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "05qn295s9kifkwr71a2rgf7y80kyrc9qli9giazav1i5dbfnk5sy";
  };
  buildPhase = "CC=gcc ./build.sh";
  installPhase = "mkdir $out && mkdir $out/bin && mv RNMC $out/bin";

}


