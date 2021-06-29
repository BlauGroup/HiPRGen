with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "e96d2861c0ee1bb210a36fdb08f3fddfb261d1b0";
  buildInputs = [gsl];
  nativeBuildInputs = [meson ninja clang sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "0vd1whx6lz0mla0b49c3mn5gfz9l3962878j73avl7ld8565vli8";
  };
  mesonBuildType = "debug";
}
