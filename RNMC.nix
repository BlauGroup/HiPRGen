with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "752d53b87ebb5668e960c924dd1a282e5d54fcf2";
  buildInputs = [gsl];
  nativeBuildInputs = [meson ninja clang sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "1cf2izrzyvwmlaj5lbhkil8hg2arzbyggjvx5fn2cr62in7y5l73";
  };
  mesonBuildType = "debug";
}
