# SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
# SPDX-License-Identifier: CC0-1.0
let
  pkgs = import (builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/3566ab7246670a43abd2ffa913cc62dad9cdf7d5.tar.gz";
    sha256 = "1wwqlc4fp00swnas5p3ipzhmgz3iyk84apg2vfglr5r4mfm3m980";
  }) { };
in
pkgs.mkShell {
  strictDeps = true;
  nativeBuildInputs = builtins.attrValues {
    inherit (pkgs)
      cargo
      clippy
      nil
      nixfmt-rfc-style
      python3
      reuse
      rust-analyzer
      rustc
      rustfmt
      singular
      ;
    inherit (pkgs.python3Packages) python-lsp-server;
  };
}
