FROM nixos/nix

RUN nix-channel --update
nix-env --install python3Packages.jupyterlab