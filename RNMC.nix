with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "4b51b8819470d028d7d43713c9609f54862f0a11";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "sha256-yn8RClfTivrwb+rCvZoEo2FVgobrouBo4Pl2wPZvaqc=";
  };
  buildPhase = ''CC="g++ -O3" ./build.sh'';
  installPhase = "mkdir $out && mkdir $out/bin && mv ./build/GMC $out/bin";

}
