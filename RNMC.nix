with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "7db64bb79daf1c1236f07d9847faf1cd3a87547b";
  buildInputs = [gsl gcc sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "06aal86lbr2mka85dhy7j838g1v38sn03hjy1g173npwi942yzmr";
  };
  buildPhase = ''CC="g++ -O3" ./build.sh'';
  installPhase = "mkdir $out && mkdir $out/bin && mv ./build/GMC $out/bin";

}
