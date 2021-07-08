with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "cd2f73c05818508347d3ab5ba2afc18510ff2c84";
  buildInputs = [gsl];
  nativeBuildInputs = [meson ninja clang sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "0xy7azk6widznzlipcz8rn55slkszyl8x17fwf45g0pxfvrgvmnv";
  };
  mesonBuildType = "debug";
}

