FROM nixos/nix

RUN nix-channel --update
RUN nix-env --install python3Packages.jupyterlab