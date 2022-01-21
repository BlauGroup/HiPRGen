FROM nixos/nix

RUN nix-channel --update
RUN echo "with import <nixpkgs> {}; [ python3 python3Packages.jupyterlab ]" > /tmp/jupyter.nix
RUN nix-env -if /tmp/jupyter.nix