with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "d643ec6651f325dda13fe45207ed4286f91134a0";
  buildInputs = [gsl];
  nativeBuildInputs = [meson ninja gcc sqlite];
  src = fetchurl {
    url = "https://github.com/danielbarter/RNMC_native/archive/${hash}.tar.gz";
    sha256 = "03zfcchng1p4d2xrdn3dykhy1lnsg0mzp2pd6dj25rdy759z6c7n";
  };
  mesonBuildType = "debug";
}
