# SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
# SPDX-License-Identifier: CC0-1.0
let
  pkgs = import (builtins.fetchTarball {
    url = "https://releases.nixos.org/nixos/unstable/nixos-25.11pre805967.62b852f6c674/nixexprs.tar.xz";
    sha256 = "sha256-0VFAX85s2QcFpcNVy7J4yn7j8Cn2IZADd1wsbT+1TYs=";
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
