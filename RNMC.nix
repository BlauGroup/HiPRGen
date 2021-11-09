with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "e1851452407555fc1385eea27a55337164b46813";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "0bk8bz0hx0w6gpw29iq8h5kxpg1y5qipvsylyi2b7jbwv20zfgmv";
  };
  buildPhase = ''CC="g++ -O3" ./build.sh'';
  installPhase = "mkdir $out && mkdir $out/bin && mv ./build/GMC $out/bin";

}
