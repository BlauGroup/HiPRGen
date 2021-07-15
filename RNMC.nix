with (import <nixpkgs> {});

stdenv.mkDerivation rec {
  name = "RNMC";
  hash = "5d9672167e439f04ff3b67aa267aba146b31f86b";
  buildInputs = [gsl];
  nativeBuildInputs = [meson ninja clang sqlite];
  src = fetchurl {
    url = "https://github.com/BlauGroup/RNMC/archive/${hash}.tar.gz";
    sha256 = "04x16z4fmx1651z1g9qf39ab9mf443jhhr1qsjdfn3fj0j2nn8gq";
  };
  mesonBuildType = "debug";
}

